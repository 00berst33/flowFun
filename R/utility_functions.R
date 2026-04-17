#' startProject
#'
#' Set working directory to a previous analysis, or create a new one.
#' ??? anything else to instantiate here ???
#'
#' @param dir_name Absolute or relative filepath giving the name of the directory
#' that you would like your analysis to be saved in. By default, a folder named
#' "Cytometry_Analysis" is created in the current working directory.
#'
#' @export
startProject <- function(dir_name = "Cytometry_Analysis") {
  if (!dir.exists(dir_name)) {
    dir.create(dir_name)
    setwd(dir_name)
  } else {
    setwd(dir_name)
  }

  # save sample info or filepath
  # save compensation matrix or filepath
  # set working directory accordingly
}


#' addMetaToTable
#'
#' @param table Table containing expression data for all samples
#' @param sample_dt Table containing sample info
#' @param join_col Column on which to join the two tables described above
#'
#'
#'
#' @return The data.table provided as input, with columns for sample information
#' added.
#' @export
addMetaToTable <- function(table, sample_dt, join_col) {
  if (join_col %in% colnames(table) && join_col %in% colnames(sample_dt)) {
    new_table <- dplyr::left_join(table, sample_dt, by = join_col)
    return(new_table)
    }
  else {
    stop(paste0("The given column name `join_col` = ", join_col, " is not present in both tables."))
  }
}


# put within GatingSet to table or something
# change name
#' @keywords internal
#' @export
getPrettyColNamesFromGatingSet <- function(input) {
  # Get cytoset if input is a GatingSet
  if (is(input, "GatingSet")) {
    input <- flowWorkspace::gs_pop_get_data(input)
  }
  # Get first element of cytoset
  table <- input[[1]] %>%
    flowCore::parameters() %>%
    Biobase::pData()

  # Get channels and markers in cytoframe
  channels <- table$name
  markers <- table$desc

  # Get new column names, in form of `channel <marker>`
  new_colnames <- sapply(seq_along(channels), function(i) {
    colname <- ifelse(is.na(markers[[i]]),
                      paste0(channels[[i]], " <", channels[[i]], ">"),
                      paste0(channels[[i]], " <", markers[[i]], ">"))
  })

  # Rename channels of each cytoframe as specified above
  # res <- lapply(seq_along(input), function (i) {
  #   lapply(seq_along(channels), function (j) {
  #     flowWorkspace::cf_rename_channel(input[[i]], channels[[j]], new_colnames[[j]])
  #     })
  #   })

  return(new_colnames)
}

# changes columns of compensation matrix to channels, markers, or both (pretty)
# and checks order of columns
#' @keywords internal
#' @export
prepareCompensationMatrix <- function(matrix, gs, pattern = " <.*") {
  # Make sure all columns are numeric
  matrix <- matrix %>%
    dplyr::select(dplyr::where(is.numeric))

  ### Make sure columns are in correct order
  # Get row indices where 1 occurs for each column
  idx <- sapply(seq_along(colnames(matrix)), function(i) {
    col <- matrix[[i]]
    idx <- match(1, col)
    return(idx)
  })

  # Assign names of vector to compensation matrix colnames
  names(idx) <- colnames(matrix)
  # Sort indices numerically
  idx <- sort(idx)

  # Reorder matrix columns according to sorted vector names
  sorted_mat <- matrix[names(idx)]

  # recall compensation matrix column names need to match those in GatingSet

  # Get column names that will be used in compensation matrix
  # if (extend_cols) {
  #   pretty <- getPrettyColNamesFromGatingSet(gs)
  #   gs_cols <- pretty
  # } else {
  #   gs_cols <- colnames(gs1[[1]])
  # }

  gs_cols <- flowWorkspace::colnames(gs[[1]])

  # Match compensation matrix column names to those in GatingSet
  match_idx <- match(colnames(sorted_mat), sub(pattern, "", gs_cols))

  if (any(is.na(match_idx))) {
    warning("Column names of provided GatingSet and compensation matrix couldn't be matched. Check `pattern` argument.")
  } else {
    result <- gs_cols[match_idx]
    colnames(sorted_mat) <- result
  }

  # Return prepared matrix
  return(sorted_mat)
}

#' gatingSetToTable
#'
#' @param gs Preprocessed GatingSet to convert
#' @param population A string giving the name of the gated population to use for
#' table creation
#'
#' @return A data.table, created from the given GatingSet
#' @export
gatingSetToTable <- function(gs, population = "root") {
  # Get cytoset with desired gated population
  cs <- flowWorkspace::gs_pop_get_data(gs, population)
  # Get sample names from GatingSet
  sn <- flowWorkspace::sampleNames(gs)

  # Get list of expression matrices for each cytoframe
  gs_list <- flowWorkspace::lapply(cs, function(gh) {flowCore::exprs(gh) %>% data.table::as.data.table()})
  # Set names of list items, so that we may create an id column when combining tables
  names(gs_list) <- sn


  # check for column name mismatches here? or GatingSet may already take care of this
  # !!! instead of requiring tidytable, maybe move the dependency to Suggests, and load it
  #   if the user has it installed

  # add metadata here using pData()


  # Bind all tables into one, by row
  table <- tidytable::bind_rows(gs_list, .id = ".id")
  # Create cell id column
  table <- table %>%
    tidytable::mutate(cell_id = seq(1, nrow(table)),
                      .after = 1)

  # Make pretty table columns
  pretty_cols <- getPrettyColNamesFromGatingSet(gs)
  match_idx <- match(colnames(table), sub(" <.*", "", pretty_cols))
  colnames(table)[!is.na(match_idx)] <- pretty_cols[stats::na.exclude(match_idx)]

  return(table)
}

#' tableToFlowSet
#'
#' Convert a table of single-cell data
#' !!! from here, add function to save .fcs files
#'
#' @param table The data table to convert to a flowSet
#' @param id_col The column containing sample IDs for the experiment
#'
#' @return A flowSet
#' @export
tableToFlowSet <- function(table, id_col = .id) {
  id_col <- rlang::enquo(id_col)

  # Get (file)names of all samples in the table
  sample_names <- table %>%
    dplyr::pull(!!id_col) %>%
    unique()

  # Make a flowFrame for each sample and put it in a list
  fs <- lapply(sample_names, function(i) {
    table %>%
      dplyr::filter(!!id_col == i) %>%
      #select(where(is.double)) %>%
      dplyr::mutate(!!id_col := match(!!id_col, sample_names)) %>% # edit channel/marker/key .csv here?
      dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>%
      flowCore::flowFrame()})
  names(fs) <- sample_names

  # Convert data type to flowSet
  fs <- methods::as(fs, "flowSet")

  return(fs)
}

flowSetToTable <- function(fs) {
  table <- fs %>%
    dplyr::ungroup() %>%
    flowCore::exprs()
}

#' flowSOMToTable
#'
#' Convert FlowSOM object to a tidytable
#'
#' @param fsom
#'
#' @return A tidytable
#' @export
flowSOMToTable <- function(fsom) {
  # Get metacluster and cluster labels
  meta_labels <- FlowSOM::GetMetaclusters(fsom)
  clust_labels <- factor(FlowSOM::GetClusters(fsom))

  # Get data.frame from FlowSOM object and make column names pretty
  dt <- fsom$data
  colnames(dt) <- fsom$prettyColnames

  # Append labels as columns to the data.table
  dt <- dt %>%
    tidytable::as_tidytable() %>%
    tidytable::mutate(.id = rep(                     # add column for sample names
                        names(fsom$metaData),
                        times = sapply(fsom$metaData, function(x) diff(x) + 1)),
                      Cluster = clust_labels,
                      Metacluster = meta_labels,
                      .keep = "all")

  # Add "clustered" attribute
  clustered <- fsom$info$parameters$colsToUse

  if (!is.null(clustered)) { # If user specified which columns to cluster on
    if (is.language(clustered)) { # if input was a variable, evaluate it
      clustered <- eval(clustered)
    }
    if (is.character(clustered)) { # if input was character, match it to prettyColNames
      match_idx <- match(clustered, names(fsom$prettyColnames))
      clustered <- fsom$prettyColnames[match_idx]
      names(clustered) <- NULL
    }
    if (is.numeric(clustered)) { # if input was numeric, use it to get appropriate column names
      clustered <- fsom$prettyColnames[clustered]
      names(clustered) <- NULL
    }
    attr(dt, "clustered") <- clustered
  } else { # If user did not specify which columns to cluster on, i.e. all columns were used
    clustered <- fsom$prettyColnames
    names(clustered) <- NULL
    attr(dt, "clustered") <- clustered
  }

  # Add "codes" attribute
  # codes <- fsom$map$codes
  # attr(dt, "codes") <- codes

  # Get needed info from FlowSOM, add as attribute
  clustering <- list("colsToUse" = fsom$info$parameters$colsToUse,
                     "clustered" = clustered,
                     "xdim" = fsom$info$parameters$xdim,
                     "ydim" = fsom$info$parameters$ydim,
                     "clustering" = fsom$metaclustering,
                     "codes" = fsom$map$codes)
  attr(dt, "clustering") <- clustering

  return(dt)
}

#' getChannelMarkerPairs
#'
#' ...
#'
#' @param ff A flowFrame or filepath to a .fcs file
#' @param save_res Boolean; should the results be saved to the current working directory
#'
#' @return A two column table specifying channel/marker pairs
#' @export
getChannelMarkerPairs <- function(ff, save_res = TRUE) {
  # If a filepath was given
  if (is.character(ff)) {
    ff <- flowCore::read.FCS(ff)
  }

  tab <- Biobase::pData(flowCore::parameters(ff)) %>%
    # select columns with marker and channel names
    dplyr::select(name, desc) %>%
    # remove channels that have no corresponding marker
    dplyr::filter(!is.na(desc))

  # Save resulting table if desired
  if (save_res) {
    utils::write.csv(tab, "channel_marker_key.csv", row.names = FALSE)
  }

  return(tab)
}

#' NEED TO HANDLE MISSING PARAMETERS FROM TRANSFORMATION
#'
#' @keywords internal
#' @export
transformTable <- function(input, transformation = NULL,
                           transform_type = c("logicle", "arcsinh", "other", "none"),
                           find_inverse = FALSE) {
  if (find_inverse) {
    # Get transformation
    temp <- attr(input, "transformation")

    if (!is.null(temp)) { # If there is a transformation attached to input
      transformation <- temp
    }
    #  <- methods::slot(transformation, "transformationId")

    # Reverse transformation if necessary
    if (transform_type == "logicle") { # make switch()?
      transformation <- flowCore::inverseLogicleTransform(transformation)
    } else if (transform_type == "arcsinh") {
      stop()
    } else if (transform_type == "other") {
      stop()
    } else if (transform_type == "none") {
      stop()
    }
  }

  # Extract list of transformations from transformList object
  transforms <- methods::slot(transformation, "transforms")

  # Apply transformation
  transformed_input <- input %>%
    tidytable::as_tidytable() %>%
    tidytable::mutate(tidytable::across(tidytable::all_of(names(transforms)),
                                        .fns = function(x) transforms[[as.character(substitute(x))]]@f(x)))

  return(transformed_input)
}

#' saveFlowSOMEmbedding
#'
#' Save essential parts of FlowSOM
#'
#' @param fsom The FlowSOM object whose mapping you would like to save
#' @param file A filepath or filename to save the resulting object under
#'
#' @export
saveFlowSOMEmbedding <- function(fsom, file) {
  fsom_minimal <- list(
    codes = fsom$map$codes,
    colsUsed = fsom$map$colsUsed,
    xdim = fsom$map$xdim,
    ydim = fsom$map$ydim,
    metaclustering = fsom$metaclustering
  )
  saveRDS(fsom_minimal, file)
}


#' loadFlowSOMEmbedding
#'
#' Load saved FlowSOM and apply to new data
#'
#' @param file A file path to the FlowSOM object whose clusters you would like `newData`
#' to map to
#' @param newData A flowFrame, flowSet, matrix with column names, or a list of
#' paths for files or directories containing the data to be clustered
#'
#' @return A FlowSOM object
#' @export
loadFlowSOMEmbedding <- function(file, newData) {
  fsom_minimal <- readRDS(file)

  # Rebuild a FlowSOM-like object shell
  fsom_shell <- FlowSOM::FlowSOM(
    input = newData,
    colsToUse = fsom_minimal$colsUsed,
    xdim = fsom_minimal$xdim,
    ydim = fsom_minimal$ydim
  )
  fsom_shell$codes <- fsom_minimal$codes

  # Project onto saved codes
  fsom_projected <- FlowSOM::NewData(fsom_shell, input = newData)

  # Attach meta-clustering back
  #fsom_projected$metaclustering <- fsom_minimal$metaclustering

  return(fsom_projected)
}


#' getSampleMetaclusterMFIs
#'
#' Get a table of MFIs, where each row is a sample and each column is a metacluster.
#'
#' @param input description
#' @param col The channel to find sample-metacluster MFIs for.
#' @param sample_df If \code{input} is a FlowSOM object, a data frame from
#' [prepareSampleInfo()].
#' @param meta_to_use A `character` list of the names of metaclusters/populations to calculate MFIs for.
#' Default is all.
#'
#' @return A table, where rows are samples and columns are metaclusters.
#' @export
getSampleMetaclusterMFIs <- function(input, col, sample_df, meta_to_use = NULL) {
  result <- UseMethod("getSampleMetaclusterMFIs")
  return(result)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.FlowSOM <- function(input, col, sample_df, meta_to_use = NULL) {
  Metacluster <- File <- NULL

  var <- rlang::enquo(col)

  # Get metacluster labels
  meta_labels <- FlowSOM::GetMetaclusters(input)

  # Append labels as columns to the data.table, and summarise
  table <- tidytable::as_tidytable(input$data) %>%
    tidytable::mutate(Metacluster = meta_labels,
                      .keep = "all") %>%
    tidytable::group_by(File) %>%
    tidytable::pivot_wider(names_from = Metacluster, values_from = !!var, id_cols = File, values_fn = median) %>%
    tidytable::filter(File %in% rownames(sample_df)) %>%
    as.data.frame()

  # change file number to sample and assign as rownames here
  rownames(table) <- sample_df[match(table$File, rownames(sample_df)), "File.Name"]

  table <- table %>%
    dplyr::select(-"File")

  if(!is.null(meta_to_use)) {
    table <- table %>%
      dplyr::select(meta_to_use)
  }

  return(table)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.data.table <- function(input, col, sample_df, meta_to_use = NULL) {
  var <- rlang::ensym(col)

  # get table where rows are samples, columns are metaclusters, and values are MFIs
  fitc_mfis <- input %>%
    dplyr::select(.id, !!var, Metacluster) %>%
    dplyr::group_by(.id, Metacluster) %>%
    dplyr::summarise(mfi = median(.data[[rlang::as_string(var)]]), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Metacluster, values_from = mfi) %>%
    dplyr::mutate(.id = basename(.id))

  fitc_mfis <- fitc_mfis %>%
    dplyr::select(!.id) %>%
    data.frame(row.names = fitc_mfis$.id, check.names = FALSE)

  if(!is.null(meta_to_use)) {
    fitc_mfis <- fitc_mfis %>%
      dplyr::select(meta_to_use)
  }

  return(fitc_mfis)
}


gs_getSampleMetaclusterMFIs <- function(gs, col, sample_df, subpopulations, inverse = FALSE) {
  var <- rlang::ensym(col)

  mfis <- flowWorkspace::gs_pop_get_stats(gs,
                                          nodes = subpopulations,
                                          type = pop.MFI,
                                          inverse = inverse) %>%
    dplyr::select(sample, pop, !!var) %>%
    tidyr::pivot_wider(names_from = pop, values_from = !!var)

  # Ensure columns are in order corresponding to design matrix
  sample_idx <- match(sample_df[, 2], mfis$sample)
  mfis <- mfis[sample_idx,]
  return(mfis)
}


#' flagMarkerNames
#'
#' Checks for typos/discrepancies between a flowSet or cytoset's marker names
#' between samples
#'
#' @param input A flowSet or cytoset
#'
#' @return A list. For each sample, contains a matrix of channel/marker pairs
#' that were inconsistent, i.e. not found in all samples
#' @export
flagMarkerNames <- function(input) {
  # Get channel/marker pairs for each sample as a list
  cm_list <- lapply(seq_len(length(input)), function(i) {cf <- input[[i]] %>% flowCore::parameters() %>%
    Biobase::pData() %>%
    dplyr::select(name, desc) %>%
    dplyr::mutate(sample_num = i)
  rownames(cf) <- NULL
  return(cf)})

  # Bind data.frames in list by row
  cm <- dplyr::bind_rows(cm_list)

  # Make table counting each time a marker name appears
  discp <- cm %>%
    dplyr::group_by(desc) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(flag = ifelse((count > num_samples | count <= 2) & !is.na(desc), 1, 0)) # add column flagging suspicious marker names

  # Get flagged marker names only
  flagged <- discp %>%
    dplyr::filter(flag == 1) %>%
    dplyr::pull(desc)

  # For each sample, subset channel/marker pairs to flagged only
  discp_list <- lapply(cm_list, function(cm) {res <- cm %>% dplyr::filter(desc %in% flagged)
                                              if (nrow(res) == 0) {res <- NULL}
                                              return(res)})

  return(discp_list)
}


#' flagChannelNames
#'
#' Checks for typos/discrepancies between a flowSet or cytoset's channel names
#' between samples
#'
#' @param input A flowSet or cytoset
#'
#' @return A list. For each sample, contains a matrix of channel/marker pairs
#' that were inconsistent, i.e. not found in all samples
#' @export
flagChannelNames <- function(input) {
  # Get channel/marker pairs for each sample as a list
  cm_list <- lapply(seq_len(length(input)), function(i) {cf <- input[[i]] %>% flowCore::parameters() %>%
    Biobase::pData() %>%
    dplyr::select(name, desc) %>%
    dplyr::mutate(sample_num = i)
  rownames(cf) <- NULL
  return(cf)})

  # Bind data.frames in list by row
  cm <- dplyr::bind_rows(cm_list)

  # Make table counting each time a channel name appears
  discp <- cm %>%
    dplyr::group_by(name) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::mutate(flag = ifelse(count > num_samples | count <= 2, 1, 0)) # add column flagging suspicious channel names

  # Get flagged channel names only
  flagged <- discp %>%
    dplyr::filter(flag == 1) %>%
    dplyr::pull(name)
  print(flagged)

  # For each sample, subset channel/marker pairs to flagged only
  discp_list <- lapply(cm_list, function(cm) {cm <- dplyr::filter(name %in% flagged)})


  discp_list <- lapply(cm_list, function(cm) {res <- cm %>% dplyr::filter(name %in% flagged)
                                              if (nrow(res) == 0) {res <- NULL}
                                              return(res)})

  return(discp_list)
}


#' makeBoolean
#'
#' @keywords internal
#' @return A boolean vector where `TRUE` means a cell belongs to the population
#' @export
makeBoolean <- function(input, indices, keep_indices = FALSE) {
  # if gatinghierarchy
  # if flowFrame
  # if cytoframe

  mat <- flowCore::exprs(input)

  row_num <- nrow(mat)

  if (keep_indices) {
    boolean_vec <- rep(FALSE, row_num)
    boolean_vec[indices] <- TRUE
  } else {
    boolean_vec <- rep(TRUE, row_num)
    boolean_vec[indices] <- FALSE
  }

  return(boolean_vec)
}

#' getClusterIndicesBySample
#'
#' @keywords internal
#' @return A list of tibbles, where each element corresponds to a sample and
#' contains the indices of cells belonging to each cluster
#' @export
getClusterIndicesBySample <- function(table) { # .id and Metacluster column assumed
  # Get list of tibbles, where each element is a sample
  sample_tibs <- table %>%
    dplyr::mutate(.id = factor(.id, levels = unique(.id))) %>%
    dplyr::group_by(.id) %>% # check column name
    dplyr::group_split() %>%
    stats::setNames(unique(table$.id))
    #stats::setNames(unique(ifelse(sapply(table$.id, is.character), basename(table$.id), table$.id)))

  # For each sample, group by metacluster and get indices of their rows (relevant to subsetted table)
  idx_tables <- lapply(sample_tibs, function(sample)
  {
    sample %>%
      dplyr::group_by(Metacluster) %>%
      dplyr::group_data()
  })

  return(idx_tables)
}

#' addClustersToGatingSet
#'
#' Uses cluster labels to add corresponding gates to the GatingSet
#'
#' @param table A data.table or tibble containing expression data for all samples.
#' Must have a column named `.id` indicating sample ID, and a column
#' `Metacluster` with cluster labels.
#' @param gs The `GatingSet` to add gates to.
#' @param parent_gate `character` indicating which population is the parent of
#' the clusters.
#'
#' This function should be used after a satisfactory clustering has been obtained
#' using the human-in-the-loop approach outlined by the vignette. The resulting
#' table should contain expression data for each sample, with each cell assigned
#' a cluster label. These labels will then be used to construct a gate for each
#' cluster.
#'
#' @export
addClustersToGatingSet <- function(table, gs, parent_gate) {
  idx_tables <- getClusterIndicesBySample(table)

  # check if gates already exist?

  # Subset cytoset to parent population
  parent_cs <- flowWorkspace::gs_pop_get_data(gs, y = parent_gate)

  # Get metacluster names
  meta_names <- idx_tables[[1]] %>%
    dplyr::pull(Metacluster) %>%
    unique()

  # Initialize list
  meta_idx <- vector("list", length(meta_names))

  # Get list of indices in metacluster for each sample
  for (i in seq_along(meta_names)) {
    # Get list of indices in metacluster for each sample
    idx_list <- sapply(idx_tables, function(table)
    {
      table %>%
        dplyr::filter(Metacluster == meta_names[i]) %>%
        dplyr::pull(.rows)
    })

    # Get list of boolean vectors for current metacluster
    bool_list <- lapply(seq_len(length(idx_list)), function(j)
    {
      makeBoolean(parent_cs[[names(idx_tables)[j]]], unlist(idx_list[[j]]), keep_indices = TRUE) # can probably remove unlist()
    }) # note keep_indices=TRUE b/c indices give cells to keep, not remove

    # Name list elements
    names(bool_list) <- names(idx_tables)

    # Assign result to list
    meta_idx[[i]] <- bool_list
  }
  # Name new list with metacluster names
  names(meta_idx) <- meta_names

  ## Add gates to GatingSet
  # !!! do you need to check if the gate already exists?
  lapply(seq_along(meta_idx), function(i) {flowWorkspace::gs_pop_add(gs, meta_idx[[i]], name = names(meta_idx)[i], parent = parent_gate)})

  ## Add keywords to GatingSet
  # Get final metaclustering
  metaclustering <- table %>%
    dplyr::select(Cluster, Metacluster) %>%
    dplyr::group_by(Cluster) %>%
    unique() %>%
    dplyr::arrange(Cluster) %>% # sort rows by ascending cluster number
    dplyr::pull(Metacluster)

  # Get metadata on clustering from data.table attributes
  clustering <- attr(table, "clustering")
  # Edit "clustering" element
  clustering[["clustering"]] <- metaclustering

  # Get first GatingHierarchy
  gh <- gs[[1]]

  # Check if any clusterings have already been added to GatingHierarchy
  # gh_clusterings <- tryCatch(
  #   keyword(gh, "clusterings"),
  #   error = function(e) {
  #    list()
  #   }
  # )

  # if (is.character(gh_clusterings)) {
  #   gh_clusterings <- list(eval(str2lang(gh_clusterings)))
  # }

  # Add new clustering to list
  #gh_clusterings[[parent_gate]] <- clustering

  # Make keyword a list
  clustering <- list(clustering)

  flowWorkspace::gh_keyword_insert(gh, parent_gate, clustering)

  # # Add edited list of clusterings to GatingSet
  # if (length(gh_clusterings) == 1) { # if this is the first clustering
  #   flowWorkspace::gh_keyword_insert(gh, "clusterings", gh_clusterings)
  # } else { # if this is a secondary clustering
  #   flowWorkspace::gh_keyword_set(gh, "clusterings", gh_clusterings)
  # }
}


# instead of this, maybe add column for proposed correction to functions above?
checkMarkerNames <- function(cm) {
  # Get channel/marker pairs for each sample as a list
  cm_list <- lapply(seq_len(length(input)), function(i) {cf <- input[[i]] %>% flowCore::parameters() %>%
    Biobase::pData() %>%
    dplyr::select(name, desc) %>%
    dplyr::mutate(position = i)
  rownames(cf) <- NULL
  return(cf)})

  # Bind data.frames in list by row
  cm <- dplyr::bind_rows(cm_list)

  # Check for typos in markers
  # probably easiest to only check for markers that don't exist; i.e.
  # CD14 may be a typo for CD16, but also is a real marker
  num_samples <- length(input)

  discp <- cm %>%
    dplyr::group_by(desc) %>%
    dplyr::summarise(count = dplyr::n()) %>%
    dplyr::filter(count < num_samples)

  likely_typo <- discp %>%
    dplyr::filter(count <= 2) %>%
    dplyr::pull(desc)
  print(likely_typo)

  likely_match <- discp %>%
    dplyr::filter(count > 2) %>%
    dplyr::pull(desc)

  distances <- adist(likely_typo, likely_match)
  colnames(distances) <- likely_match

  # for each row, find similar marker names that are most likely to be a match
  closest_match <- lapply(seq_len(nrow(distances)), function(i) {
    row <- distances[i,]
    min_dist <- min(row)
    closest_match <- names(row[row==min_dist])
    return(closest_match)})
  names(closest_match) <- likely_typo

  return(closest_match)
}


###
# get vector of marker names
# make boolean vector based on whether or not they were used for clustering
  # know if they were used for clustering based on vector of clustered markers given
#
# markers <- markernames(cs)
#
# old_pdata <- pData(cs)
# num_samples <- nrow(old_pdata)
# # for marker in markers
# # initialize vector rep(FALSE, nrow(old_pdata))
# # if marker is in cols_to_cluster (and maybe add check for sample) then change value to TRUE
# clustered <- lapply(markers, function(marker) {
#   if (marker %in% cols_to_cluster) { # To check for samples, could loop over them instead of using rep()
#     init <- rep(TRUE, num_samples) # need to check if cols_to_cluster is character vector
#   } else {
#     init <- rep(FALSE, num_samples)
#   }
#   df <- data.frame(init)
#   colnames(df) <- marker
#   return(df)
#   })
#
# clustered <- dplyr::bind_cols(result)
# result <- cbind(old_pdata, clustered)


# input is GatingSet
addMarkersToMetaData <- function(input, cols_to_cluster) {

  markers <- markernames(input)

  old_pdata <- pData(input)
  num_samples <- nrow(old_pdata)

  clustered <- lapply(markers, function(marker) {
    if (marker %in% cols_to_cluster) { # To check for samples, could loop over them instead of using rep()
      init <- rep(TRUE, num_samples) # need to check if cols_to_cluster is character vector
    } else {
      init <- rep(FALSE, num_samples)
    }
    df <- data.frame(init)
    colnames(df) <- marker
    return(df)
  })

  clustered <- dplyr::bind_cols(clustered)
  result <- cbind(old_pdata, clustered)

  # edit actual pData
  # pData(input) <- result
  return(result)

  # no need to return anything in final
}


#' addMetadataToGatingSet
#'
#'
#' @param gs A `GatingSet`
#' @param sample_df A `data.frame`, where each row corresponds to a sample and
#' each column corresponds to metadata about the samples, such as experimental group,
#' or the filename of a corresponding control sample. Must contain a column for
#' filename, called `File.Name`.
#'
#' @export
addMetadataToGatingSet <- function(gs, sample_df) {
  # add control files names to metadata
  cs <- gs_pop_get_data(gs)

  # Get metadata from cytoset
  cs_meta <- cs %>%
    pData()

  # Match filenames in sample information and metadata tables
  idx <- match(sample_df$File.Name, rownames(cs_meta))
  # Order sample info, then get all metadata excluding filenames
  new_metadata <- sample_df[idx, -"File.Name"]

  # Add to metadata table
  cs_meta <- cbind(cs_meta, new_metadata)

  # Set pData
  flowCore::pData(cs) <- cs_meta
}
