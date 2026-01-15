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
  id_col <- enquo(id_col)

  # Get (file)names of all samples in the table
  sample_names <- table %>%
    pull(!!id_col) %>%
    unique()

  # Make a flowFrame for each sample and put it in a list
  fs <- lapply(sample_names, function(i) {
    table %>%
      filter(!!id_col == i) %>%
      #select(where(is.double)) %>%
      mutate(!!id_col := match(!!id_col, sample_names)) %>% # edit channel/marker/key .csv here?
      select(where(is.numeric)) %>%
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
    select(name, desc) %>%
    # remove channels that have no corresponding marker
    filter(!is.na(desc))

  # Save resulting table if desired
  if (save_res) {
    write.csv(tab, "channel_marker_key.csv", row.names = FALSE)
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
#' @param csv_name Optional. If you would like the resulting table to be saved
#' to a .csv file, define the desired filename here.
#'
#' @return A table, where rows are samples and columns are metaclusters.
#' @export
getSampleMetaclusterMFIs <- function(input, col, csv_name = NULL) {
  result <- UseMethod("getSampleMetaclusterMFIs")
  return(result)
}

#' @keywords internal
#' @export
getSampleMetaclusterMFIs.FlowSOM <- function(input, col, sample_df, csv_name = NULL) {
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
    write.csv(table, csv_name)
  }

  table <- table %>%
    dplyr::select(-"File")

  return(table)
}

# returns list of tables for data table
getSampleMetaclusterMFIsNew <- function(input, sample_df, cols_to_use,
                                        linear = FALSE, meta_names = NULL,
                                        ctrl_input = NULL, subsetted_meta = NULL) {
  Metacluster <- .id <- NULL

  # Set a seed for reproducibility
  set.seed(42)

  # Get metacluster names if none were given, and filter out those not of interest
  if (!is.null(meta_names)) {
    input <- input %>%
      tidytable::filter(Metacluster %in% meta_names)
  } else {
    meta_names <- input %>%
      tidytable::pull(Metacluster) %>%
      unique()
  }

  # Get all file names
  filenames <- input %>%
    tidytable::pull(.id) %>%
    unique() %>%
    basename()

  # Determine samples excluded from analysis, if any
  removed_samples <- which(!(filenames %in% sample_df[, 2]))

  # Master data frame that will contain data about all metaclusters.
  df_full <- data.frame()

  # Get the channel corresponding to the current control
  channels <- NULL
  fcs_cols <- attr(input, "cols_from_fcs")

  for (i in seq_along(fcs_cols)) {
    col <- fcs_cols[i]

    if (col %in% cols_to_use) {
      channels <- c(channels, col)
    } else {
      if (!is.null(attr(input[[col]], "marker"))) {
        if (!is.na(attr(input[[col]], "marker")))
          if (attr(input[[col]], "marker") %in% cols_to_use) {
            channels <- c(channels, col)
          }
      }
    }

  }
  if (is.null(channels)) {
    stop(paste0("No corresponding channels were found for ", cols_to_use,
                " check parameter `cols_to_use` for typos."))
  }

  # Transform data back to linear scale if desired, or if applying controls
  if (linear & !is.null(attr(input, "transformation"))) {
    input <- transformTable(input, find_inverse = TRUE)
  }

  # Get sample/metacluster MFIs for each channel in main input
  input <- input %>%
    tidytable::mutate(.id = basename(.id)) %>%
    tidytable::rename(file = .id) %>%
    tidytable::summarise(tidytable::across(.cols = fcs_cols,
                                           .fns = stats::median,
                                           .drop = "keep"),
                         .by = c(file, Metacluster))

  ### DEALING WITH CTRL
  # assume `ctrl_input` is clustered data table
  # double check conditionals here and when re-assigning `input`

  # for (col in cols_to_use) {
  #   # If channel of interest has a corresponding control given
  #   if (!is.null(ctrl_input)) { # check if it has a control
  #     # Transform to linear scale
  #     ctrl_input <- transformTable(ctrl_input, find_inverse = TRUE)
  #
  #     # Get sample/metacluster MFIs for current marker/channel of interest
  #     ctrl_input <- ctrl_input %>%
  #       tidytable::filter(Metacluster %in% meta_names) %>%
  #       tidytable::mutate(.id = basename(.id)) %>%
  #       tidytable::rename(file = .id) %>%
  #       tidytable::summarise(tidytable::across(.cols = channel,
  #                                              .fns = stats::median,
  #                                              .drop = "keep"),
  #                            .by = c(file, Metacluster)) %>%
  #       tidytable::pivot_wider(names_from = Metacluster, values_from = channel)
  #   } # then subtract these values from corresponding sample/metacluster pairs
  # }

  # Get sample/metacluster MFIs for marker/channels of interest
  tables <- lapply(channels, function(channel) {
    input %>%
      tidytable::select(c(file, Metacluster, channel)) %>%
      tidytable::pivot_wider(names_from = Metacluster, values_from = channel)})

  # Bind tables together, and add a row for marker/channel
  tables <- tidytable::bind_rows(tables, .id = TRUE) %>%
    tidytable::mutate(channel = channels[.id],
                      .keep = "unused",
                      .after = 2)

  ##### bits from original function
  # Checks for samples with too few cells
  # missing_samples <- which(!(rownames(sample_df) %in% input))

  return(tables)
}

doDEAnalysisNew <- function(input, sample_df, design, contrasts, cols_to_use,
                            meta_names = NULL, ctrl_input = NULL, subsetted_meta = NULL,
                            save_csv = FALSE, dir_tables = NULL) {
  value <- channel <- NULL

  # Create empty data frames for each comparison that will be made
  pval_dfs <- lapply(1:ncol(contrasts), function(i) {data.frame()})

  # Create expression matrix for testing, if it wasn't given
  if (!all(c("file", "channel") %in% colnames(input))) {
    input <- getSampleMetaclusterMFIsNew(input, sample_df, cols_to_use, meta_names = meta_names)
  }
  expr_matrix <- input %>%
    tidytable::pivot_longer(cols = tidytable::where(is.numeric),
                            names_to = "Metacluster",
                            values_to = "value") %>%
    tidytable::pivot_wider(names_from = file,
                           values_from = value)

  channel_order <- expr_matrix %>%
    tidytable::pull(channel) %>%
    unique()

  # Turn input into a matrix
  expr_matrix <- expr_matrix %>%
    tidytable::select(-1) %>%
    as.matrix(rownames = "Metacluster")

  # Create linear models.
  # NOTE: weights and logic questionable !!!
  lm_model <- limma::lmFit(object = expr_matrix,
                           design = design)

  # Perform statistical tests.
  contrasts_fit <- limma::contrasts.fit(lm_model, contrasts)
  limma_ebayes <- limma::eBayes(contrasts_fit, trend = TRUE)

  # If a directory to write tables in wasn't given
  if (is.null(dir_tables)) {
    dir_tables <- dir_limma()
  }

  # Create tables containing results of our statistical tests, and add them
  # to the data frame corresponding to the relevant comparison.
  for (i in seq_len(ncol(contrasts))) {
    table <- limma::topTable(limma_ebayes,
                             coef = colnames(contrasts)[i],
                             sort.by = "p",
                             number = Inf,
                             adjust.method = "BH")

    # Add column for marker ID
    if (nrow(table) != 0) {
      table$marker <- rep(cols_to_use[1], nrow(table))

      if (length(cols_to_use) > 1) { # if more than one marker is being tested
        for (j in 1:nrow(table)) {
          curr <- rownames(table)[j]
          ix <- ceiling(as.integer(curr)/length(meta_names))
          table$marker[j] <- channel_order[ix]
        }
        rownames(table) <- NULL
      }

      # Append table of results to returned list
      temp_ind <- nrow(pval_dfs[[i]]) + 1
      pval_dfs[[i]] <- rbind(pval_dfs[[i]], table)
    }

    # Save table, if desired
    if (save_csv) {
      utils::write.csv(pval_dfs[[i]], file = file.path(dir_limma(),
                                                       dir_tables,
                                                       paste0(colnames(contrasts)[i],
                                                              "_", dir_tables, "_limma.csv")))
    }
  }

  return(pval_dfs)
}
