#' getMarker
#'
#' Helper for getting a flowFrame channel's corresponding marker.
#'
#' @param ff A flowFrame.
#' @param channel The channel of interest.
#'
#' @return A marker name.
#'
#' @export
getMarker <- function(ff, channel) {
  nm <- methods::as(Biobase::pData(flowCore::parameters(ff))[
    which(Biobase::pData(flowCore::parameters(ff))[, "name"] == channel), "desc"], "character")
  if (is.na(nm)) {
    return(channel)
  } else {
    return(nm)
  }
}

#' getChannel
#'
#' Helper for getting a flowFrame marker's corresponding channel.
#'
#' @inheritParams getMarker
#' @param marker The marker of interest.
#'
#' @return A channel name.
#'
#' @export
getChannel <- function(ff, marker) {
  return(methods::as(Biobase::pData(flowCore::parameters(ff))[
    which(Biobase::pData(flowCore::parameters(ff))[, "desc"] == marker), "name"], "character"))
}

#' plotBeforeAfter
#'
#' Plot events that were removed between steps.
#'
#' @param ff1 The "before" flowFrame.
#' @param ff2 The "after" flowFrame.
#' @param channel1 The channel to plot on the x-axis.
#' @param channel2 The channel to plot on the y-axis.
#' @param ncells The number of cells to include in the plot.
#'
#' @return A plot where cells that were removed are highlighted in red, and cells
#' that were kept are highlighted in blue.
#'
#' @export
plotBeforeAfter <- function(ff1, ff2, channel1, channel2, ncells) {
  x <- y <- NULL

  df <- data.frame(x = flowCore::exprs(ff1)[, channel1],
                   y = flowCore::exprs(ff1)[, channel2])
  i <- sample(nrow(df), ncells)
  if (!"Original_ID" %in% colnames(flowCore::exprs(ff1))) {
    ff1@exprs <- cbind(ff1@exprs,
                       Original_ID = seq_len(nrow(ff1@exprs)))
  }
  plot <- ggplot2::ggplot(df[i,], ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(size = 0.5,
                        color = ifelse(flowCore::exprs(ff1)[i,"Original_ID"] %in%
                                         flowCore::exprs(ff2)[,"Original_ID"], 'blue', 'red')) +
    ggplot2::xlab(getMarker(ff1, channel1)) +
    ggplot2::ylab(getMarker(ff1, channel2)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  ggpubr::ggarrange(plot)
}

#' getTableFromFCS
#'
#' Generate a data.table to use for analysis from existing .fcs files.
#'
#' @param input List of filenames with relative file paths, or a directory name.
#'
#' @return A data.table containing .fcs file data.
#'
#' @export
getTableFromFCS <- function(input) {
  # Prepare input
  if (!is.list(input) && dir.exists(input)) {
    files <- list.files(input, full.names = TRUE)
  } else {
    files <- input
  }

  if (length(files) == 1) {
    # Read in file
    ff <- flowCore::read.FCS(files, truncate_max_range = FALSE)

    data <- flowCore::exprs(ff)

    prepr_table <- tidytable::data.table(.id = rep(files, nrow(data)), cell_id = seq(1, nrow(data)), data)

  } else if (length(files) > 1) { # check
    # Initialize final data table
    prepr_tables <- lapply(1:length(files), function(i) {tidytable::data.table()})

    # First cell ID
    current_id <- 1

    # Iterate over all files
    for (i in 1:length(files)) {
      file <- files[[i]]

      # Read in file
      ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)

      # Make cell ID column
      vec <- seq(current_id, current_id + nrow(flowCore::exprs(ff)) - 1)

      # Add data table to list
      dt <- tidytable::data.table(cell_id = vec, flowCore::exprs(ff))
      prepr_tables[[i]] <- rbind(prepr_tables[[i]], dt)

      # Set first cell ID for next file
      current_id <- nrow(flowCore::exprs(ff)) + 1
    }
    # Concatenate all data tables into one, with column for sample ID
    data.table::setattr(prepr_tables, 'names', files)
    prepr_table <- tidytable::bind_rows(prepr_tables, .id = TRUE) %>%
      tidytable::select(dplyr::where(~ !any(is.na(.)))) # remove columns that don't appear in all tables
  }

  return(prepr_table)
}


#' doPreprocessing
#'
#' Preprocess all files in the given directory.
#' !!! WIP edit gates and handling of low-quality samples, rename markers if necessary
#'
#' @param dir_raw The directory containing raw data files to be processed together.
#' @param ld_channel Name of the channel corresponding to the marker for live/dead
#' cells. This should appear the same as it does in the .fcs files.
#' @param compensation A compensation matrix, or a file path to one.
#' @param transformation The transformation to apply to the data.
#' @param debris_gate A gate to gate out debris. If \code{NULL} (default), a
#' gate is automatically selected.
#' @param live_gate A gate to gate out dead cells. If \code{NULL} (default), a
#' gate is automatically selected.
#' @param nmad Parameter to determine strictness of doublet removal. See
#' [PeacoQC::RemoveDoublets()] for details.
#' @param pctg_live Minimum percentage of live cells a sample should have.
#' @param pctg_qc Minimum percentage of cells a sample should have after removing
#' low-quality events during quality control.
#' @param pdf_name Name of the PDF containing diagnostic plots for preprocessing
#' results. Default is \code{"preprocessing_results.pdf"}.
#' @param save_fcs Boolean, should the preprocessed data be saved as .fcs files.
#' Default is \code{FALSE}.
#'
#' @details
#' The steps taken by this function, in order, are:
#'   1. Removing margin events (i.e. values that fall outside the range the
#'      flow cytometer should be able to pick up).
#'   2. Removing doublets.
#'   3. Removing debris.
#'   4. Compensating the data.
#'   5. Performing a log-icle transformation on each channel.
#'   6. Removing dead cells.
#'   7. Quality control, to check and correct for occurrences like clogs, speed
#'      changes, etc. and flag any files that may be too low-quality to include
#'      in the analysis.
#'
#' The script will output a few files to check for any mistakes. The first is
#' a PDF highlighting the cells that were removed in each .fcs file. The others
#' are plots for any files that were flagged during quality control, where each
#' channel is plotted against time and removed events are marked.
#'
#' @return A data.table containing preprocessed data.
#'
#' @export
doPreprocessing <- function(dir_raw, ld_channel, compensation, transformation,
                            debris_gate = NULL, live_gate = NULL, nmad = 4,
                            pctg_live = 0.6, pctg_qc = 0.8,
                            pdf_name = "preprocessing_results.pdf", save_fcs = FALSE) {
  Plot <- df <- NULL

  # Get all files
  raw_files <- list.files(dir_raw, full.names = TRUE)

  # Get all measured channels in flowFrame, excluding FSC and SSC
  all_channels <- as.vector(Biobase::pData(flowCore::parameters(ff))$name[which(
    !is.na(Biobase::pData(flowCore::parameters(ff))$desc))])

  prepr_tables <- lapply(1:length(raw_files), function(i) {tidytable::data.table()})

  # Preprocess all given files and generate PDF of preprocessing results
  grDevices::pdf(file.path("Preprocessing Results", pdf_name))
  for (file in raw_files) {
    print(file)
    ff <- flowCore::read.FCS(file.path(dir_raw, file), truncate_max_range = FALSE)

    # RENAME HERE

    # Remove margin events and doublets
    ff_m <- PeacoQC::RemoveMargins(ff, c("FSC-A", all_channels))
    ff_d <- PeacoQC::RemoveDoublets(ff_m, nmad = nmad)

    # Create plot for results of doublet removal
    p1 <- plotBeforeAfter(ff_m, ff_d, "FSC-A", "FSC-H", 5000)

    # Apply debris gate if given, otherwise draw automatically
    if (!is.null(debris_gate)) {
      ff_g <- ff_d[flowCore::filter(ff_d, debris_gate)@subSet, ]
    } else {
      stop() # temp
    }

    # Create plot for results of debris removal
    p2 <- plotBeforeAfter(ff_d, ff_g, "FSC-A", "SSC-A", 5000)

    # Compensate and transform data
    ff_c <- flowCore::compensate(ff_g, compensation)
    ff_t <- flowCore::transform(ff_c, transformation)

    # Apply live/dead gate if given, otherwise draw automatically
    if (!is.null(live_gate)) {
      ff_l <- ff_t[flowCore::filter(ff_t, live_gate)@subSet, ]
    } else {
      stop() # temp
    }

    # Create plot for results of live/dead cell removal
    p3 <- plotBeforeAfter(ff_t, ff_l, ld_channel, "FSC-A", 5000)

    # Arrange above plots on single page
    gridExtra::grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = file)

    # Perform quality control via flowCut package
    fc <- suppressWarnings(flowCut::flowCut(f = ff_l,
                                            FileID = file,
                                            Plot = Plot,
                                            Directory = file.path("Preprocessing Results",
                                                                  paste0(dir_raw, " QC")),
                                            Verbose = TRUE))
    ff_fc <- fc$frame

    # Print a warning for any files with less than 60% live cells.
    if (nrow(ff_l)/nrow(ff_g) < pctg_live) {
      warning(paste0(file, " has only ", round(nrow(ff_l)/nrow(ff_g)*100, 2),
                     "% live cells."))
    }
    # Print a warning for any files that had more than 20% of its events removed by QC.
    if (nrow(ff_fc)/nrow(ff_l) < pctg_qc) {
      warning(paste0(file, " had ", round(nrow(ff_fc)/nrow(ff_l), 2),
                     " of its events removed by flowCut."))
    }

    # Bind to greater data table
    dt <- tidytable::data.table(flowCore::exprs(ff_fc))
    prepr_tables[[file]] <- rbind(prepr_tables[[file]], dt)

    if (save_fcs) {
      flowCore::write.FCS(ff_fc, file = paste0("Preprocessed ", dir_raw, file))
    }

  }
  grDevices::dev.off()

  # Concatenate all data tables into one, with column for sample ID
  prepr_table <- data.table::rbindlist(prepr_tables, use.names = FALSE, idcol = TRUE)
  return(prepr_table)
}
