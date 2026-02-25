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
  table <- tidytable::bind_rows(gs_list, .id = "File")
  # Create cell id column
  table <- table %>%
    tidytable::mutate(cell_id = seq(1, nrow(table)),
                      .after = 1)

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

  # Append labels as columns to the data.table
  dt <- fsom$data
  dt <- dt %>%
    tidytable::as_tidytable() %>%
    tidytable::mutate(Cluster = clust_labels,
                      Metacluster = meta_labels,
                      .keep = "all")

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
#' @param csv_name Optional. If you would like the resulting table to be saved
#' to a .csv file, define the desired filename here.
#'
#' @return A table, where rows are samples and columns are metaclusters.
#' @export
getSampleMetaclusterMFIs <- function(input, col, sample_df, csv_name = NULL) {
  result <- UseMethod("getSampleMetaclusterMFIs")
  return(result)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.FlowSOM <- function(input, col, sample_df, csv_name = NULL) {
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

  if(!is.null(csv_name)) {
    utils::write.csv(table, csv_name)
  }

  table <- table %>%
    dplyr::select(-"File")

  return(table)
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
  discp_list <- lapply(cm_list, function(cm) {cm <- dplyr::filter(desc %in% flagged)})

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

  # For each sample, subset channel/marker pairs to flagged only
  discp_list <- lapply(cm_list, function(cm) {cm <- dplyr::filter(name %in% flagged)})

  return(discp_list)
}



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
