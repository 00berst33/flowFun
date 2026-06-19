#' Function to match a character vector whose elements are either markers or
#' channels to a parent vector whose elements contain both markers and
#' the corresponding channel. i.e. match "CD3" to "CD3 <Alexa Fluor 700-A>".
#'
#' @keywords internal
#' @export
getFullNames <- function(substr, target_names, keep_names = FALSE) {
  # Helper to escape regex metacharacters
  escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }

  # Get new names of input
  full_names <- sapply(substr, function(x) {
    x_esc <- escape_regex(x)

    # token-style matching
    pattern <- paste0(
      "(?<![[:alnum:]])",
      x_esc,
      "(?![[:alnum:]])"
    )
    grep(pattern, target_names, value = TRUE, perl = TRUE)
  })


  if (!keep_names) {
    # Remove names given by grep()
    full_names <- unname(full_names)
  }
  full_names <- unlist(full_names)

  return(full_names)
}


#' addMetaToTable
#'
#' @param table Table containing expression data for all samples
#' @param sample_dt Table containing sample info
#' @param join_col Column on which to join the two tables described above. Or,
#' to join on different variables between tables, an expression (e.g. `a == b`).
#'
#'
#' @return The data.table provided as input, with columns for sample information
#' added.
#' @export
addMetaToTable <- function(table, sample_dt, join_col) {
  var <- rlang::enquo(join_col)

  new_table <- dplyr::left_join(table, sample_dt, by = dplyr::join_by(!!var))
  return(new_table)
}


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


#' prepareCompensationMatrix
#'
#' Changes columns of compensation matrix to channels, markers, or both (pretty)
#' and checks order of columns
#'
#' @param matrix The compensation matrix to prepare
#' @param gs The `GatingSet` whose parameter names the compensation matrix's
#' columns should be matched to
#'
#' @export
prepareCompensationMatrix <- function(matrix, gs) {
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

  # compensation matrix column names need to match those in GatingSet
  # Get GatingSet column names
  gs_cols <- flowWorkspace::colnames(gs[[1]])

  # Match compensation matrix column names to those in GatingSet
  # Assume gs_cols is regular channel names
  new_colnames <- getFullNames(gs_cols, colnames(sorted_mat), keep_names = TRUE) %>%
    names()

  if (any(is.na(new_colnames))) {
    warning("Column names of provided GatingSet and compensation matrix couldn't be matched. Check that column names of the GatingSet are substrings of the compensation matrix's column names.")
  } else {
    colnames(sorted_mat) <- new_colnames
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

  # pull marker <channel> into two different vector with regexpr if there are <>

  # Get (file)names of all samples in the table
  sample_names <- table %>%
    dplyr::pull(!!id_col) %>%
    unique()

  # Make a flowFrame for each sample and put it in a list
  fs <- lapply(sample_names, function(i) {
    ff <- table %>%
      dplyr::filter(!!id_col == i) %>%
      dplyr::select(-!!id_col) %>%
      dplyr::select(dplyr::where(is.numeric)) %>%
      as.matrix() %>%
      flowCore::flowFrame()
    flowCore::keyword(ff)[["FIL"]] <- i
    return(ff)})
  names(fs) <- sample_names

  # Convert data type to flowSet
  fs <- methods::as(fs, "flowSet")

  return(fs)
}


#' @keywords internal
#' @export
flowSetToTable <- function(fs) {
  table <- lapply(seq_along(fs), function(i) {
    df <- flowCore::exprs(fs[[i]]) %>%
      as.data.frame(make.names = FALSE) %>%
      dplyr::mutate(.id = flowCore::sampleNames(fs)[[i]],
                    .before = 1)
  }) %>%
    dplyr::bind_rows() %>%
    tidytable::as_tidytable()

  return(table)
}


#' @keywords internal
#' @export
flowSetMakePrettyNames <- function(fs) {
  channels <- flowCore::colnames(fs)
  ff <- fs[[1]]

  pretty_names <- sapply(seq_along(channels), function(i) {
    marker <- getMarker(ff, channels[i])
    pretty_name <- paste0(marker, " <", channels[i], ">")
  })

  return(pretty_names)
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
inverseTransformTable <- function(input, transformation) {

  # Reverse transformation
  transforms <- lapply(transformation, `[[`, "inverse")

  if (is.null(transforms)) {
    stop("No inverse found for `transformation`.")
  }

  # If column names don't match
  if (!any(names(transforms) %in% colnames(input))) {
    # Extract only channel names from input
    input_channels <- sub(".*<(.*)>.*", "\\1", colnames(input))
    if (all(names(transforms) %in% input_channels)) {
      match_idx <- match(names(transforms), input_channels)
      # Set transform names to pretty column names
      names(transforms) <- colnames(input)[match_idx]
    } else {
      stop("One or more names in the transformation and `input` don't match.")
    }
  }

  # Get channel names in transformation
  channels <- names(transforms)

  # Transform relevant channels in data.table
  for (i in seq_along(channels)) {

    channel <- channels[[i]]
    trans_fun <- transforms[[i]]

    input <- input %>%
      tidytable::as_tidytable() %>%
      tidytable::mutate(
        !!channel := tidytable::map_dbl(.data[[channel]], trans_fun)
      )
  }

  return(input)
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
#' @param input A `data.table` or `GatingSet`
#' @param col The channel to find sample-metacluster MFIs for.
#' @param populations A `character` list of the names of metaclusters/populations to calculate MFIs for.
#' @param transformation A \code{\link[flowCore:transformList]{transformList}} specifying
#' the transformation applied to the data, if any. Required when `input` is a `data.table`
#' and `inverse = TRUE`
#' @param inverse Boolean, whether data should be back-transformed before calculating MFIs.
#' Only valid when `input` is a `GatingSet` containing a transformation.
#'
#' @return A table, where rows are samples and columns are metaclusters.
#' @export
getSampleMetaclusterMFIs <- function(input, col, populations = NULL,
                                     transformation = NULL, inverse = FALSE) {
  result <- UseMethod("getSampleMetaclusterMFIs")
  return(result)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.GatingSet <- function(input, col, populations = NULL,
                                               transformation = NULL, inverse = FALSE) {
  var <- rlang::ensym(col)

  # if (is.null(populations)) {
  #   stop("Must specify `populations` parameter when `input` is a GatingSet!")
  # }

  # Get MFIs for populations of interest
  mfis <- flowWorkspace::gs_pop_get_stats(input, nodes = populations, type = pop.MFI, inverse = inverse)

  # If column is a channel, find associated marker
  if (!(col %in% colnames(mfis))) {
    cf <- flowWorkspace::gh_pop_get_data(input[[1]])
    pdata <- flowWorkspace::pData(flowCore::parameters(cf))
    idx <- match(col, pdata$name)
    col <- pdata[idx, "desc"]
  }

  # If column is a channel, find associated marker
  if (!(col %in% colnames(mfis))) {
    cf <- flowWorkspace::gh_pop_get_data(input[[1]])
    pdata <- flowWorkspace::pData(flowCore::parameters(cf))
    idx <- match(col, pdata$name)
    col <- pdata[idx, "desc"]
  }

  # Select marker/channel of interest and pivot data to wide format
  mfis <- mfis %>%
    dplyr::select(sample, pop, !!col) %>%
    tidyr::pivot_wider(names_from = pop, values_from = !!col)

  # Set rownames
  rn <- mfis$sample
  mfis <- mfis %>%
    as.data.frame(make.names = FALSE) %>%
    dplyr::select(-sample)
  rownames(mfis) <- rn

  return(mfis)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.data.table <- function(input, col, transformation,
                                                inverse = FALSE, populations = NULL) {
  # Back transform data, if necessary
  if (inverse & is.null(transformation)) {
    stop("Must specify `transformation` when `inverse = TRUE` for `data.tables`!")
  } else if (inverse & !is.null(transformation)) {
    input <- inverseTransformTable(input, transformation)
  }

  # If input is only channel name
  if (is.character(col) & !(col %in% colnames(input))) {
    marker_cols <- sub(".*<(.*)>.*", "\\1", colnames(input))
    match_idx <- match(col, marker_cols)
    # Set col to pretty column names
    col <- colnames(input)[match_idx]
  }

  var <- rlang::ensym(col)

  # get table where rows are samples, columns are metaclusters, and values are MFIs
  mfis <- input %>%
    dplyr::select(.id, !!var, Metacluster) %>%
    dplyr::group_by(.id, Metacluster) %>%
    dplyr::summarise(mfi = median(.data[[rlang::as_string(var)]]), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Metacluster, values_from = mfi) %>%
    dplyr::mutate(.id = basename(.id))

  mfis <- mfis %>%
    dplyr::select(!.id) %>%
    data.frame(row.names = mfis$.id, check.names = FALSE)

  if (!is.null(populations)) {
    mfis <- mfis %>%
      dplyr::select(populations)
  }

  return(mfis)
}


#' makeBoolean
#'
#' Helper for addClustersToGatingSet
#'
#' @keywords internal
#' @return A boolean vector where `TRUE` means a cell belongs to the population
#' @export
makeBoolean <- function(input, indices, keep_indices = FALSE) {
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
#' Helper for addClustersToGatingSet
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
#' @param fsom_file Optional, a character giving a path to an .rds file for a FlowSOM
#' object. Default is `NULL`, which checks for a filename in `attr(table, "fsom_filename")`
#'
#' This function should be used after a satisfactory clustering has been obtained
#' using the human-in-the-loop approach outlined by the vignette. The resulting
#' table should contain expression data for each sample, with each cell assigned
#' a cluster label. These labels will then be used to construct a gate for each
#' cluster.
#'
#' @export
addClustersToGatingSet <- function(table, gs, parent_gate, fsom_file = NULL) {
  # Get cluster row indices for each sample
  # Returns a list of tables, one for each sample
  idx_tables <- getClusterIndicesBySample(table)

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
      makeBoolean(parent_cs[[names(idx_tables)[j]]], unlist(idx_list[[j]]), keep_indices = TRUE)
    }) # note keep_indices=TRUE b/c indices give cells to keep, not remove

    # Name list elements
    names(bool_list) <- names(idx_tables)

    # Assign result to list
    meta_idx[[i]] <- bool_list
  }
  # Name new list with metacluster names
  names(meta_idx) <- meta_names

  ## Add gates to GatingSet
  lapply(seq_along(meta_idx), function(i) {
    tryCatch({flowWorkspace::gs_pop_add(gs, meta_idx[[i]], name = names(meta_idx)[i], parent = parent_gate)},
             error = function(e) { # If gate already exists, remove it and redraw it
               message(paste0("The gate ", names(meta_idx)[[i]], " already exists. It will be deleted and redrawn."))
               flowWorkspace::gs_pop_remove(gs, node = names(meta_idx)[i])
               flowWorkspace::gs_pop_add(gs, meta_idx[[i]], name = names(meta_idx)[i], parent = parent_gate)
             }) # recompute?
    })

  ## Add keywords to GatingSet
  # Get final metaclustering
  metaclustering <- table %>%
    dplyr::select(Cluster, Metacluster) %>%
    dplyr::group_by(Cluster) %>%
    unique() %>%
    dplyr::arrange(Cluster) %>% # sort rows by ascending cluster number
    dplyr::pull(Metacluster)

  if (is.null(fsom_file)) { # if fsom_file is NULL, check if there is a filename in table attributes
    fsom_file <- attr(table, "fsom_filename")
  }
  if (!is.null(fsom_file)) {
    flowWorkspace::pData(gs)[, parent_gate] <- fsom_file

    # Edit FlowSOM file
    fsom <- readRDS(fsom_file)
    fsom$metaclustering <- metaclustering
    saveRDS(fsom, file = fsom_file)
  }
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


#' addMetadataToGatingSet
#'
#'
#' @param gs A `GatingSet`
#' @param sample_df A `data.frame`, where each row corresponds to a sample and
#' each column corresponds to metadata about the samples, such as experimental group,
#' or the filename of a corresponding control sample. Must contain a column for
#' filename, called `filename`.
#'
#' @export
addMetadataToGatingSet <- function(gs, sample_df) {
  # add control files names to metadata
  cs <- gs_pop_get_data(gs)

  # Get metadata from cytoset
  cs_meta <- cs %>%
    pData()

  # Match filenames in sample information and metadata tables
  idx <- match(sample_df$filename, rownames(cs_meta))
  # Order sample info, then get all metadata excluding filenames
  new_metadata <- sample_df[idx, ] %>%
    dplyr::select(-filename)

  # Add to metadata table
  cs_meta <- cbind(cs_meta, new_metadata)

  # Set pData
  flowCore::pData(cs) <- cs_meta
}


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
