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

    data <- tidytable::as_tidytable(flowCore::exprs(ff))

    prepr_table <- data %>%
      tidytable::mutate(.id = rep(files, nrow(data)),
                        cell_id = seq(1, nrow(data)),
                        .before = 1)
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

  # Get all channel names from .fcs files
  fcs_channels <- as.vector(Biobase::pData(flowCore::parameters(ff))$name)
  fcs_col_nums <- which(fcs_channels %in% c(colnames(prepr_table)))

  # Add attribute to table specifying which columns came from .fcs files
  data.table::setattr(prepr_table, "cols_from_fcs", fcs_channels[fcs_col_nums])
  # Add attributes to keep track of which marker each channel/column corresponds to
  for (i in fcs_col_nums) {
    col <- fcs_channels[i]
    data.table::setattr(prepr_table[[col]], "marker", as.vector(Biobase::pData(flowCore::parameters(ff))$desc)[i])
  }

  return(prepr_table)
}

# implementation of postgres
# give user option to re-select gate manually
# break down into smaller functions
# let user preview results of gates


#' doPreprocessing
#'
#' Preprocess all files in the given directory.
#' !rename markers if necessary
#' !save setting used somehow, .json, .csv, etc.
#'
#' @param input A directory, list of filenames, or data.table.
#' @param ld_channel Name of the channel corresponding to the marker for live/dead
#' cells. This should appear the same as it does in the .fcs files.
#' @param compensation Optional. A compensation matrix, or a file path to one.
#' @param transformation The transformation to apply to the data. If \code{NULL}
#' (default), a log-icle transformation is applied to the data.
#' @param transformation_type The type of transformation to be applied to the data.
#' If \code{transformation} is \code{NULL}, the transformation parameters will be
#' selected automatically. See details for more.
#' @param debris_gate A gate to gate out debris, defined as a
#' \code{\link[flowCore:filter-class]{filter}}. If \code{NULL} (default), a
#' gate is automatically selected with \code{\link[flowDensity:flowDensity-methods]{flowDensity()}}.
#' @param live_gate A gate to gate out dead cells. If \code{NULL} (default), a
#' gate is automatically selected.
#' @param nmad Parameter to determine strictness of doublet removal. See
#' [PeacoQC::RemoveDoublets()] for details.
#' @param pctg_live Minimum proportion of live cells a sample should have remaining
#' after gating out dead cells. The sample will be excluded if this number isn't met.
#' @param pctg_qc Minimum proportion of cells a sample should have remaining
#' after low-quality events are removed during quality control. The sample will be
#' excluded if this number isn't met.
#' @param pdf_name Name of the PDF containing diagnostic plots for preprocessing
#' results. Default is \code{"preprocessing_results.pdf"}.
#' @param flowcut_dir The name of the directory containing QC results. The plots
#' for any flagged files will be saved here. Default is \code{"flowCut"}.
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
#' a PDF highlighting the cells that were removed in each .fcs file, whose name is
#' specified by \code{pdf_name}. The others are plots for any files that were
#' flagged during quality control, written in the directory specified by
#' \code{flowCut}. Each channel is plotted against time and removed events are marked.
#'
#' The user may define their own transformation to apply to the data, if
#' desired. This may be done by creating an object of type
#' \code{\link[flowCore:transformList-class]{transformList}}, containing a transformation
#' for each channel of interest, and passing it to the function parameter
#' \code{transformation}. The parameter \code{transformation_type} should then
#' be selected accordingly.
#'
#' \code{transformation_type} specifies the type of transformation to apply to the data:
#' \code{"logicle"}, \code{"arcsinh"}, \code{"other"}, or \code{"none"}. Regardless of
#' whether or not \code{transformation} is \code{NULL}, it is necessary to specify
#' this parameter so that the original values may be retrieved later. So, if
#' a custom hyperbolic sine transformation was given, this parameter should be set
#' to \code{"arcsinh"}, and if a custom log transformation was given, it
#' should be set to \code{"other"}, etc.
#'
#' If NO custom transformation was given, and \code{"logicle"} was chosen, then
#' the function will automatically determine appropriate parameters for the
#' transformation and apply it. If \code{"arcsinh"} was chosen, a hyperbolic sine
#' transformation with cofactor 5 will be selected and applied. Finally, if
#' \code{"none"} is chosen, the data will have no transformation applied. This is
#' generally not advised, as it will make most data visualizations used
#' by this workflow difficult or impossible to interpret.
#'
#' For flow data, a log-icle transformation is highly recommended, as it represents
#' low signals and compensated data better than a typical logarithmic scaling,
#' and results in more interpretable data visualizations. For CyTOF data, a
#' hyperbolic sine transformation with a cofactor of 5 is standard, and also recommended.
#'
#' The user may also define their own debris and live/dead gates to apply to their samples
#' via the parameters \code{debris_gate} and \code{live_gate}. The objects passed to
#' these parameters should be of class \code{\link[flowCore:filter-class]{filter}}.
#' If these gates are not specified by the user, gates will be automatically drawn using
#' \code{\link[flowDensity:flowDensity-methods]{flowDensity()}}.
#'
#' Note that appropriate compensation and transformation of the data is crucial
#' for an accurate live/dead gate. The step of drawing said gate may be
#' skipped by leaving the parameter \code{ld_channel} unspecified.
#'
#' Samples will be removed from the returned table if they do not meet the proportion of
#' live and quality cells specified by \code{pctg_live} and \code{pctg_qc}, as a sample
#' with too few viable events is likely not worth including in the analysis. However,
#' if it is preferred for no samples to be removed, setting both parameters to \code{0}
#' will result in all samples being included in the returned table, regardless of quality.
#'
#' @return A data.table containing preprocessed data.
#'
#' @seealso [PeacoQC::RemoveMargins()], [PeacoQC::RemoveDoublets()], [flowCore::compensate()],
#' [flowCore::transform()], [flowCore::estimateLogicle()], [flowCore::logicleTransform()],
#' [flowCore::arcsinhTransform()], [flowCore::biexponentialTransform()],
#' [flowDensity::deGate()], [flowCut::flowCut()]
#'
#' @export
doPreprocessing <- function(input, ld_channel = NULL, compensation = NULL, transformation = NULL,
                            transformation_type = c("logicle", "arcsinh", "other", "none"),
                            debris_gate = NULL, live_gate = NULL, nmad = 4,
                            pctg_live = 0.6, pctg_qc = 0.8,
                            pdf_name = "preprocessing_results.pdf", flowcut_dir = "flowCut",
                            save_fcs = FALSE) {
  Plot <- df <- .id <- NULL

  # Prepare data
  if (!methods::is(input, "data.table")) {
    if (file.exists(input) | dir.exists(input) | is.list(input)) {
      input <- getTableFromFCS(input)
    } else {
      print("`input` should be either an existing directory, a list of filenames, or a table generated by `getTableFromFCS()`.")
    }
  }

  # Get all files
  raw_files <- input %>%
    tidytable::pull(.id) %>%
    unique()

  # Get all measured channels in flowFrame, excluding FSC and SSC
  fcs_cols <- attr(input, "cols_from_fcs") # which have marker attribute that isn't NA
  all_channels <- c()
  for (i in seq_along(fcs_cols)) {
    col <- fcs_cols[i]
    if(!is.na(attr(input[[col]], "marker"))) {
      all_channels <- c(all_channels, col)
    }
  }

  # Initialize list of data.tables to store results
  prepr_tables <- lapply(1:length(raw_files), function(i) {tidytable::data.table()})

  # Create directory if needed
  dir <- "Preprocessing Results"
  if (!dir.exists(dir)) {
    dir.create(dir)
  }

  # Initialize gating scheme
  gating_scheme <- list("debris_gate" = list("channels" = c("FSC-A", "SSC-A"),
                                             "filter" = "flowDensity"),
                        "live_gate" = list("channels" = c(ld_channel, "FSC-A"),
                                           "filter" = "flowDensity")) # allow for user specified gates, later

  debris_percent <- 0.85
  live_percent <- pctg_live

  # Change gate types, if needed
  if (inherits(debris_gate, "polygonGate")) {
    debris_gate <- methods::slot(debris_gate, "boundaries")
    gating_scheme[["debris_gate"]][["filter"]] <- debris_gate
  }
  if (inherits(live_gate, "polygonGate")) {
    live_gate <- methods::slot(live_gate, "boundaries")
    gating_scheme[["live_gate"]][["filter"]] <- live_gate
  }

  # Preprocess all given files and generate PDF of preprocessing results
  grDevices::pdf(file.path(dir, pdf_name))
  for (file in raw_files) {
    print(paste0("Processing ", basename(file), "..."))
    ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)

    # RENAME HERE

    # Remove margin events and doublets
    ff_m <- PeacoQC::RemoveMargins(ff, c("FSC-A", all_channels))
    ff_d <- PeacoQC::RemoveDoublets(ff_m, nmad = nmad)

    # Create plot for results of doublet removal
    p1 <- plotBeforeAfter(ff_m, ff_d, "FSC-A", "FSC-H", 5000)

    # Apply debris gate if given, otherwise draw automatically
    if (!is.null(debris_gate)) {
      ff_g <- flowDensity::flowDensity(obj = ff_d,
                                       channels = c("FSC-A", "SSC-A"),
                                       position = c(TRUE, NA),
                                       filter = debris_gate)
      ff_g <- flowDensity::getflowFrame(ff_g)
    } else {
      ff_g <- flowDensity::flowDensity(obj = ff_d,
                                       channels = c("FSC-A", "SSC-A"),
                                       position = c(TRUE, NA),
                                       percentile = c(debris_percent, NA)) # base on last gate if it exists
      debris_percent <- methods::slot(ff_g, "proportion")/100
      ff_g <- flowDensity::getflowFrame(ff_g)
    }

    # Create plot for results of debris removal
    p2 <- plotBeforeAfter(ff_d, ff_g, "FSC-A", "SSC-A", 5000)

    # Apply compensation matrix if given
    if (!is.null(compensation)) {
      ff_c <- flowCore::compensate(ff_g, spillover = compensation)
    } else {
      ff_c <- ff_g
    }

    # Apply transformation
    if (!is.null(transformation)) {
      ff_t <- flowCore::transform(ff_c, transformation)
    } else {
      transformation <- switch(transformation_type,
                               "logicle" = flowCore::estimateLogicle(ff_c, channels = all_channels),
                               "arcsinh" = {
                                 res <- flowCore::arcsinhTransform(b = 5)
                                 res <- flowCore::transformList(all_channels, res)
                                 res},
                               "none" = {
                                 res <- flowCore::linearTransform()
                                 res <- flowCore::transformList(all_channels, res)
                                 res
                               },
                               "other" = stop("Must provide a value for `transformation` when `transformation_type` is 'other'."))

      ff_t <- flowCore::transform(ff_c, transformation)
    }

    # If there is a live/dead channel
    if (!is.null(ld_channel)) {
      # Apply live/dead gate if given, otherwise draw automatically
      if (!is.null(live_gate)) {
        ff_l <- flowDensity::flowDensity(obj = ff_t,
                                         channels = c("ld_channel", "FSC-A"),
                                         position = c(FALSE, NA),
                                         filter = live_gate)
        ff_l <- flowDensity::getflowFrame(ff_g)
      } else {
        ff_l <- flowDensity::flowDensity(obj = ff_t,
                                         channels = c(ld_channel, "FSC-A"),
                                         position = c(FALSE, NA),
                                         percentile = c(live_percent, NA),
                                         twin.factor = c(0.1, NA))
        live_percent <- methods::slot(ff_l, "proportion")/100
        ff_l <- flowDensity::getflowFrame(ff_l)
      }
    }

    # Create plot for results of live/dead cell removal
    p3 <- plotBeforeAfter(ff_t, ff_l, ld_channel, "FSC-A", 5000)

    # Arrange above plots on single page
    gridExtra::grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = file)

    # Perform quality control via flowCut package
    fc <- suppressWarnings(flowCut::flowCut(f = ff_l,
                                            FileID = make.names(file),
                                            Plot = "Flagged Only",
                                            Directory = file.path("Preprocessing Results",
                                                                  "flowCut"),
                                            Verbose = TRUE))
    ff_fc <- fc$frame

    remove <- FALSE

    # Print a warning for any files with less than the given proportion of live cells.
    if (nrow(ff_l)/nrow(ff_g) < pctg_live) {
      warning(paste0(file, " has only ", round(nrow(ff_l)/nrow(ff_g)*100, 2),
                     "% live cells."))
      remove <- TRUE
    }
    # Print a warning for any files that had more than the given proportion of events removed by QC.
    if (nrow(ff_fc)/nrow(ff_l) < pctg_qc) {
      warning(paste0(file, " had ", 1 - round(nrow(ff_fc)/nrow(ff_l)*100, 2),
                     "% of its events removed by flowCut."))
      remove <- TRUE
    }

    # If it is high enough quality, bind sample to greater data table and save
    if (!remove) {
      dt <- tidytable::data.table(flowCore::exprs(ff_fc))
      prepr_tables[[file]] <- rbind(prepr_tables[[file]], dt)

      if (save_fcs) {
        flowCore::write.FCS(ff_fc, file = paste0("Preprocessed ", file))
      }
    } else {
      warning(paste0(file, " was removed due to poor quality. If desired, parameters `pctg_live`",
      " and `pctg_qc` may be edited to change the strictness of sample removal."))
    }
  }
  grDevices::dev.off()

  # Concatenate all data tables into one, with column for sample ID
  prepr_table <- tidytable::bind_rows(prepr_tables, .id = TRUE)
  prepr_table <- prepr_table %>%
    tidytable::mutate(cell_id = seq(1, nrow(prepr_table)),
                      .before = 2)

  # Add attributes
  if (is.null(compensation)) {
    compensation <- "none"
  }
  data.table::setattr(prepr_table, "compensation_matrix", compensation)

  data.table::setattr(prepr_table, "transformation", transformation)
  methods::slot(attr(prepr_table, "transformation"), "transformationId") <- transformation_type

  data.table::setattr(prepr_table, "gating_scheme", gating_scheme)

  # Get all channel names from .fcs files
  fcs_channels <- as.vector(Biobase::pData(flowCore::parameters(ff))$name)
  fcs_col_nums <- which(fcs_channels %in% colnames(prepr_table))

  # Add attribute to table specifying which columns came from .fcs files
  data.table::setattr(prepr_table, "cols_from_fcs", fcs_channels[fcs_col_nums])
  # Add attributes to keep track of which marker each channel/column corresponds to
  for (i in fcs_col_nums) {
    col <- fcs_channels[i]
    data.table::setattr(prepr_table[[col]], "marker", as.vector(Biobase::pData(flowCore::parameters(ff))$desc)[i])
  }

  return(prepr_table)
}


#' prepareGateParams
#'
#' @keywords internal
#' @param gate A gate defined by either the flowDensity or flowCore packages, or a list.
#'
#' @export
prepareGateParams <- function(gate) {
  # Set default options
  default_options <- list(position = c(TRUE, NA))

  if (inherits(gate, "parameterFilter")) {
    gate <- methods::as(gate, "polygonGate")
    boundaries <- methods::slot(gate, "boundaries")
    channels <- colnames(boundaries)

    all_options <- utils::modifyList(default_options, list("channels" = channels,
                                                           "filter" = boundaries))
  } else if (inherits(gate, "CellPopulation")) {
    filter <- methods::slot(gate, "filter")
    channels <- colnames(filter)

    all_options <- utils::modifyList(default_options, list("channels" = channels,
                                                           "filter" = filter))
  } else if (is.list(gate)) {
    # Add any additional options chosen by user, and overwrite defaults if needed
    all_options <- utils::modifyList(default_options, gate)
  } else {
    stop("The provided gate's data type is not supported.")
  }

  return(all_options)
}

#' previewPreprocessing
#'
#' Preview preprocessing settings for a single file before applying it to all data.
#'
#' @inheritParams doPreprocessing
#' @param file A filepath or flowFrame.
#'
#' @export
previewPreprocessing <- function(file, ld_channel = NULL, compensation = NULL, transformation = NULL,
                                 transformation_type = c("logicle", "arcsinh", "other", "none"),
                                 debris_gate = NULL, live_gate = NULL, nmad = 4) {
  if (file.exists(file)) {
    ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
  } else {
    ff <- file
  }

  all_channels <- as.vector(Biobase::pData(flowCore::parameters(ff))$name[which(!is.na(Biobase::pData(flowCore::parameters(ff))$desc))])

  # RENAME HERE

  # Prepare parameters for gating functions
  debris_defaults <- list("channels" = c("FSC-A", "SSC-A"),
                          "position" = c(TRUE, NA),
                          "percentile" = c(0.85, NA))
  if (is.null(debris_gate)) {
    debris_gate <- prepareGateParams(debris_defaults)
  } else if (is.list(debris_gate)) {
    debris_gate <- prepareGateParams(utils::modifyList(debris_defaults, debris_gate))
  } else {
    debris_gate <- prepareGateParams(debris_gate)
  }

  if (!is.null(ld_channel) & !is.null(compensation)) {
    live_defaults <- list("channels" = c(ld_channel, "SSC-A"),
                          "position" = c(FALSE, NA),
                          "percentile" = c(0.65, NA),
                          "twin.factor" = c(0.1, NA))
    if (is.null(live_gate)) {
      live_gate <- prepareGateParams(live_defaults)
    } else if (is.list(live_gate)) {
      live_gate <- prepareGateParams(utils::modifyList(live_defaults, live_gate))
    } else {
      live_gate <- prepareGateParams(live_gate)
    }
  }

  # Remove margin events and doublets
  ff_m <- PeacoQC::RemoveMargins(ff, c("FSC-A", "SSC-A", all_channels))
  ff_d <- PeacoQC::RemoveDoublets(ff_m, nmad = nmad)

  # Create plot for results of doublet removal
  p1 <- plotBeforeAfter(ff_m, ff_d, "FSC-A", "FSC-H", 5000)

  # Apply debris gate
  ff_g <- do.call(flowDensity::flowDensity, c(ff_d, debris_gate))
  ff_g <- flowDensity::getflowFrame(ff_g)

  # Create plot for results of debris removal
  p2 <- plotBeforeAfter(ff_d, ff_g, "FSC-A", "SSC-A", 5000)

  # Apply compensation matrix if given
  if (!is.null(compensation)) {
    ff_c <- flowCore::compensate(ff_g, spillover = compensation)
  } else {
    ff_c <- ff_g
  }

  # Apply transformation
  if (!is.null(transformation)) {
    ff_t <- flowCore::transform(ff_c, transformation)
  } else {
    transformation <- switch(transformation_type,
                             "logicle" = flowCore::estimateLogicle(ff_c, channels = all_channels),
                             "arcsinh" = {
                               res <- flowCore::arcsinhTransform(b = 5)
                               res <- flowCore::transformList(all_channels, res)
                               res},
                             "none" = {
                               res <- flowCore::linearTransform()
                               res <- flowCore::transformList(all_channels, res)
                               res
                             },
                             "other" = stop("Must provide a value for `transformation` when `transformation_type` is 'other'."))

    ff_t <- flowCore::transform(ff_c, transformation)
  }

  # Apply live/dead gate
  if (!is.null(live_gate)) {
    ff_l <- do.call(flowDensity::flowDensity, c(ff_t, live_gate))
    ff_l <- flowDensity::getflowFrame(ff_l)
  }

  # Create plot for results of live/dead cell removal
  p3 <- plotBeforeAfter(ff_t, ff_l, ld_channel, "FSC-A", 5000)

  # Arrange above plots
  gridExtra::grid.arrange(p1, p2, p3, nrow = 2, ncol = 2, top = file)

  return(ff_l)
}

