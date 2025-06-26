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
  i <- sample(nrow(df), min(nrow(df), ncells))
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

  #ggpubr::ggarrange(plot)
  return(plot)
}

#' getTableFromFCS
#'
#' Generate a data.table to use for analysis from existing .fcs files.
#'
#' @param input List of filenames with relative file paths, or a directory name.
#' @param num_cells Optional parameter for stratified sampling, the number of
#' cells to randomly sample from each file.
#'
#' @return A data.table containing .fcs file data.
#'
#' @export
getTableFromFCS <- function(input, num_cells = NULL) {
  # Prepare input
  if (all(dir.exists(input))) {
    # files <- lapply(input, list.files, full.names = TRUE)
    files <- unlist(lapply(input, list.files, full.names = TRUE), use.names = FALSE)
  } else if (all(file.exists(input))) {
    files <- input
  } else {
    stop("Invalid input.")
  }

  if (is.null(num_cells)) {
    num_cells <- Inf
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

      # Make cell ID column and add to table
      vec <- seq(current_id, current_id + nrow(flowCore::exprs(ff)) - 1)
      dt <- tidytable::data.table(cell_id = vec, flowCore::exprs(ff))

      # Set first cell ID for next file
      current_id <- nrow(flowCore::exprs(ff)) + 1

      # Sample `num_cells`
      ids <- sample.int(length(vec), size = min(num_cells, length(vec)))
      ids <- sort(ids)
      dt <- dt[ids, ]

      # Add data table to list
      prepr_tables[[i]] <- rbind(prepr_tables[[i]], dt)
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
# avoid opening several pdfs when the function keeps failing

#' doPreprocessing
#'
#' Preprocess all files in the given directory.
#' !rename markers if necessary
#'
#' @param input A directory, list of filenames, or data.table.
#' @param num_cells The number of cells to sample from each sample. Default is
#' 50000. To use as many cells as possible, this parameter may be set to \code{Inf}.
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
#' @param ld_channel Name of the channel corresponding to the marker for live/dead
#' cells. This should appear the same as it does in the .fcs files.
#' @param nmad Parameter to determine strictness of doublet removal. See
#' [PeacoQC::RemoveDoublets()] for details.
#' @param pctg_live Minimum proportion of live cells a sample should have remaining
#' after gating out dead cells. The sample will be excluded if this number isn't met.
#' @param pctg_qc Minimum proportion of cells a sample should have remaining
#' after low-quality events are removed during quality control. The sample will be
#' excluded if this number isn't met.
#' @param save_plots Boolean, should plots for gating and quality control
#' results be saved. Default is \code{TRUE}.
#' @param save_fcs Boolean, should the preprocessed data be saved as .fcs files.
#' Default is \code{FALSE}.
#' @param pdf_name Name of the PDF containing diagnostic plots for preprocessing
#' results. Default is \code{"preprocessing_results.pdf"}.
#' @param flowcut_dir The name of the directory containing QC results. The plots
#' for any flagged files will be saved here. Default is \code{"flowCut"}.
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
doPreprocessing <- function(input, num_cells = 50000, compensation = NULL, transformation = NULL,
                            transformation_type = c("logicle", "arcsinh", "other", "none"),
                            debris_gate = NULL, ld_channel = NULL, live_gate = NULL, nmad = 4,
                            pctg_live = 0.6, pctg_qc = 0.8,
                            save_plots = TRUE, save_fcs = FALSE,
                            pdf_name = "preprocessing_results.pdf"
                            #flowcut_dir = "flowCut"
                            ) {
  Plot <- df <- .id <- NULL

  # Prepare data
  # if (!methods::is(input, "data.frame")) {
  #   if (is.vector(input) | all(file.exists(input)) | all(dir.exists(input))) {
  #     input <- getTableFromFCS(input, num_cells = Inf)
  #   } else {
  #     print("`input` should be either an existing directory, a list of filenames, or a table generated by `getTableFromFCS()`.")
  #   }
  # }

  if (!methods::is(input, "flowSet")) {
      if (!methods::is(input, "data.frame")) {
        if (is.vector(input) | all(file.exists(input)) | all(dir.exists(input))) {
          input <- getTableFromFCS(input, num_cells = Inf)
        }
      }
    raw_files <- input %>%
      tidytable::pull(.id) %>%
      unique()

    # Get all measured channels in flowFrame, excluding FSC and SSC
    fcs_cols <- attr(input, "cols_from_fcs") # which have marker attribute that isn't NA
    all_channels <- c()
    for (i in seq_along(fcs_cols)) {
      col <- fcs_cols[i]
      marker <- attr(input[[col]], "marker")
      if(!is.na(marker)) {
        all_channels <- c(all_channels, col)
      }
    }
  } else if (methods::is(input, "flowSet")) {
    raw_files <- input
    all_channels <- as.vector(Biobase::pData(flowCore::parameters(input[[1]]))$name[which(!is.na(Biobase::pData(flowCore::parameters(input[[1]]))$desc))])
  } else if (!methods::is(input, "data.frame")) {
      print("`input` should be either an existing directory, a list of filenames, or a table generated by `getTableFromFCS()`.")
  }

  # Initialize list of data.tables to store results
  prepr_tables <- lapply(1:length(raw_files), function(i) {data.table::data.table()})

  # Initialize list to store plots
  grobs <- list()

  # Create directory if needed
  # dir <- "Preprocessing Results"
  # if (!dir.exists(dir)) {
  #   dir.create(dir)
  # }

  # Initialize gating scheme
  gating_scheme <- list()

  # Set default parameters for debris gate
  debris_defaults <- list("channels" = c("FSC-A", "SSC-A"),
                          "position" = c(TRUE, NA),
                          "percentile" = c(0.85, NA))

  # Add additional parameters specified by user, if necessary
  if (is.null(debris_gate)) {
    debris_gate <- debris_defaults
  } else if (is.list(debris_gate)) {
    debris_gate <- utils::modifyList(debris_defaults, debris_gate)
  } else {
    debris_gate <- prepareGateParams(debris_gate)
  }
  gating_scheme <- append(gating_scheme, list("debris_gate" = debris_gate))

  # If a compensation matrix was given
  if (!is.null(compensation)) {
    # Save to current project directory
    write.csv(compensation, "compensation_matrix.csv")
    # If there is a marker for live/dead cells
    if(!is.null(ld_channel)) {
      # Set default parameters for live/dead gate
      live_defaults <- list("channels" = c(ld_channel, "SSC-A"),
                            "position" = c(FALSE, NA),
                            "percentile" = c(0.70, NA),
                            "twin.factor" = c(0.1, NA))
      # Add additional parameters specified by user, if necessary
      if (is.null(live_gate)) {
        live_gate <- live_defaults
      } else if (is.list(live_gate)) {
        live_gate <- utils::modifyList(live_defaults, live_gate)
      } else {
        live_gate <- prepareGateParams(live_gate)
      }
      gating_scheme <- append(gating_scheme, list("live_gate" = live_gate))
    }
  }

  # Preprocess all given files and generate PDF of preprocessing results
  for (i in seq_along(raw_files)) {
    if (methods::is(raw_files, "flowSet")) {
      ff <- raw_files[[i]]
      file <- flowCore::identifier(ff)
    } else {
      file <- raw_files[i]
      print(paste0("Processing ", basename(file), "..."))
      ff <- flowCore::read.FCS(file, truncate_max_range = FALSE)
    }

    # RENAME HERE

    # Remove margin events and doublets
    ff_m <- PeacoQC::RemoveMargins(ff, c("FSC-A", all_channels))
    ff_d <- PeacoQC::RemoveDoublets(ff_m, nmad = nmad)

    # Apply debris gate
    ff_g <- do.call(flowDensity::flowDensity, c(list("obj" = ff_d), debris_gate))
    debris_percent <- methods::slot(ff_g, "proportion")/100
    ff_g <- flowDensity::getflowFrame(ff_g)

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

    if (!is.null(live_gate)) {
      ff_l <- do.call(flowDensity::flowDensity, c(list("obj" = ff_t), live_gate))
      live_percent <- methods::slot(ff_l, "proportion")/100
      ff_l <- flowDensity::getflowFrame(ff_l)
    }

    # Create plot for results of doublet removal
    p1 <- plotBeforeAfter(ff_m, ff_d, "FSC-A", "FSC-H", 5000)

    # Create plot for results of debris removal
    p2 <- plotBeforeAfter(ff_d, ff_g, "FSC-A", "SSC-A", 5000)

    # Create plot for results of live/dead cell removal
    p3 <- plotBeforeAfter(ff_t, ff_l, ld_channel, "FSC-A", 5000)

    j = 1 + 3*(i-1)
    grobs[[j]] <- ggplot2::ggplotGrob(p1)
    grobs[[j+1]] <- ggplot2::ggplotGrob(p2)
    grobs[[j+2]] <- ggplot2::ggplotGrob(p3)

    # Perform quality control via flowCut package
    fc <- flowCut::flowCut(f = ff_l,
                          FileID = make.names(file),
                          Plot = "All",
                          Directory = file.path(#"Preprocessing Results",
                                                "flowCut"),
                          Verbose = TRUE)
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
      dt <- data.table::data.table(flowCore::exprs(ff_fc))
      ids <- sample.int(nrow(dt), size = min(num_cells, nrow(dt))) # sample if needed
      dt <- dt[ids, ]
      prepr_tables[[i]] <- rbind(prepr_tables[[i]], dt)
      # print(as.data.frame(prepr_tables[[i]]))

      if (save_fcs) {
        flowCore::write.FCS(ff_fc, file = paste0("Preprocessed ", file))
      }
    } else {
      warning(paste0(file, " was removed due to poor quality. If desired, parameters `pctg_live`",
                     " and `pctg_qc` may be edited to change the strictness of sample removal."))
    }
  }

  if (methods::is(raw_files, "flowSet")) {
    raw_files <- Biobase::sampleNames(raw_files)
    names(prepr_tables) <- raw_files
  } else {
    names(prepr_tables) <- raw_files
  }

  # Save PDF of gating plots for each sample if desired
  if (save_plots) {
    names(grobs) <- basename(raw_files)
    #grDevices::pdf(file.path(dir, pdf_name))
    ml <- gridExtra::marrangeGrob(grobs,
                                  nrow = 2,
                                  ncol = 2,
                                  layout_matrix = rbind(c(1,1,2,2),
                                                        c(1,1,2,2),
                                                        c(NA,3,3,NA),
                                                        c(NA,3,3,NA)),
                                  top = quote(names(grobs)[g]))
    ggplot2::ggsave(pdf_name, plot = ml, device = grDevices::pdf)
    #grDevices::dev.off()
  }

  # Concatenate all data tables into one, with column for sample ID
  prepr_table <- tidytable::bind_rows(prepr_tables, .id = TRUE)
  #if (!("cell_id" %in% colnames(prepr_table))) {
  prepr_table <- prepr_table %>%
    tidytable::mutate(cell_id = seq(1, nrow(prepr_table)),
                      .before = 2)

  # if (!is.null(orig_ids)) {
  #   prepr_table <- prepr_table %>%
  #     tidytable::filter(which(cell_id %in% orig_ids))
  # }
  # }

  # Add attributes
  if (is.null(compensation)) {
    compensation <- "none"
  }

  # CHANGE TO SAVING OR MAKING LIST OBJECT HERE
  ###
  saveRDS(transformation, "transformation.rds")
  saveRDS(gating_scheme, "gating_scheme.rds")
  # save table in case of multiple clusterings with same preprocessing steps?
  ###
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
    # temp_col <- which(colnames(prepr_table) == col)
    # colnames(prepr_table)[temp_col] <- paste0(col, " <", attr(prepr_table[[col]], "marker"), ">")
  }
  prepr_table[]
  return(prepr_table)
}


#' doPreprocessing.GatingSet
#'
#' @param input A \code{GatingSet} created with the \code{flowWorkspace} package
#' @param compensation A numeric matrix
#' @param transformation A character string specifying which transformation to apply,
#' or a \code{transformList} defining a custom one
#' @param debris_args A string defining any additional gating arguments to use
#' for the non-debris gate; see [openCyto::gate_mindensity()]
#' @param singlet_args A string defining any additional gating arguments to use
#' for the non-debris gate; see [flowStats::gate_singlet()]
#' @param live_args A string defining any additional gating arguments to use
#' for the live cell gate; see [openCyto::gate_mindensity()]
#' @param ld_channel Optional; a string specifying the name of the channel used
#' to detect live/dead cells
#'
#' @return No return value, the input GatingSet is edited
#' @export
doPreprocessing.GatingSet <- function(input,
                                      compensation = NULL,
                                      transformation = c("logicle", "arcsinh", "linear"),
                                      debris_args = list(),
                                      singlet_args = list(),
                                      live_args = list(),
                                      ld_channel = NULL) {
  Plot <- df <- .id <- NULL # ???


  # If a compensation matrix was given
  if (!is.null(compensation)) {
    # Apply compensation matrix
    input <- flowCore::compensate(input, comp_mat)
  }

  # Apply transformation if it was given/desired
  if (!is.null(transformation)) {
    all_channels <- names(flowWorkspace::markernames(input))
    if(is.character(transformation)) {
      transformation <- switch(transformation,
                               "logicle" = flowCore::estimateLogicle(input[[1]], channels = all_channels),
                               "arcsinh" = {
                                 res <- flowCore::arcsinhTransform(b = 5)
                                 res <- flowCore::transformList(all_channels, res)
                                 res},
                               "linear" = {
                                 res <- flowCore::linearTransform()
                                 res <- flowCore::transformList(all_channels, res)
                                 res
                               })
    }
    input <- flowCore::transform(input, transformation)
  }

  # Debris
  openCyto::gs_add_gating_method(input,
                                 alias = "nonDebris",
                                 pop = "+",
                                 parent = "root",
                                 dims = "FSC-A",
                                 gating_method = "mindensity",
                                 gating_args = debris_args) # string

  # Doublets
  openCyto::gs_add_gating_method(input,
                                 alias = "singlets",
                                 pop = "+",
                                 parent = "nonDebris",
                                 dims = "FSC-A,FSC-H",
                                 gating_method = "singletGate",
                                 gating_args = singlet_args)

  # If there is a marker for live/dead cells
  if(!is.null(ld_channel)) {
    if(is.null(compensation)) { # !!! is compensation matrix essential? I don't remember
      warning("`ld_channel` specified but no compensation matrix was given.")
    }
    # Add live/dead gate to Gating Set
    openCyto::gs_add_gating_method(input,
                                   alias = "live",
                                   pop = "-",
                                   parent = "singlets",
                                   dims = paste0(ld_channel),
                                   gating_method = "mindensity",
                                   gating_args = live_args)
  }


  # Save PDF of gating plots for each sample if desired
  ### make this separate function?
  # if (save_plots) {
  #   pl <- lapply(input, function(gh) {ggcyto::ggcyto_arrange(ggcyto::autoplot(gh, bins = 250, strip.text = "gate"), nrow = 2)})
  #   names(pl) <- sampleNames(input)
  #   ml <- gridExtra::marrangeGrob(pl,          # should change to uniquely created plots rather than lapply()
  #                                 nrow = 1,
  #                                 ncol = 1,
  #                                 top = quote(names(pl)[g]))
  #   ggplot2::ggsave(pdf_name, plot = ml, device = grDevices::pdf, pointsize = 8, scale = 2)
  # }
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

  # If `gate` was drawn with flowCore
  if (inherits(gate, "parameterFilter")) {
    gate <- methods::as(gate, "polygonGate")
    boundaries <- methods::slot(gate, "boundaries")
    channels <- colnames(boundaries)

    all_options <- utils::modifyList(default_options, list("channels" = channels,
                                                           "filter" = boundaries))
  } else if (inherits(gate, "CellPopulation")) {
    filters <- methods::slot(gate, "filter")
    channels <- colnames(filters)

    all_options <- utils::modifyList(default_options, list("channels" = channels,
                                                           "filter" = filters))
  } else if (is.list(gate)) {
    # Add any additional options chosen by user, and overwrite defaults if needed
    all_options <- utils::modifyList(default_options, gate)
  }
  else {
    stop("The provided gate's data type is not supported. Please ensure that it is either a
         gate drawn with the flowCore/flowDensity packages, or a list.")
  }

  return(all_options)
}


#' previewPreprocessing
#'
#' Preview preprocessing settings for a single file before applying it to all data.
#'
#' @inheritParams doPreprocessing
#' @param input A filepath or flowFrame.
#'
#' @export
previewPreprocessing <- function(input, ld_channel = NULL, compensation = NULL, transformation = NULL,
                                 transformation_type = c("logicle", "arcsinh", "other", "none"),
                                 debris_gate = NULL, live_gate = NULL, nmad = 4) {
  .id <- NULL
  if (methods::is(input, "data.frame")) {
    input <- input %>%
      tidytable::pull(.id) %>%
      unique()
    input <- input[1]
    print("Checking first file found in table..")
  }

  # Read in input, if necessary
  if (file.exists(input)) {
    ff <- flowCore::read.FCS(input, truncate_max_range = FALSE)
  } else {
    ff <- input
  }

  all_channels <- as.vector(Biobase::pData(flowCore::parameters(ff))$name[which(!is.na(Biobase::pData(flowCore::parameters(ff))$desc))])

  # RENAME HERE

  # Set default parameters for debris gate
  debris_defaults <- list("channels" = c("FSC-A", "SSC-A"),
                          "position" = c(TRUE, NA),
                          "percentile" = c(0.85, NA))

  # Add to or override defaults with user specified parameters
  if (is.null(debris_gate)) {
    debris_gate <- debris_defaults
  } else if (is.list(debris_gate)) {
    debris_gate <- utils::modifyList(debris_defaults, debris_gate)
  } else {
    debris_gate <- prepareGateParams(debris_gate)
  }

  # If there is a marker for live/dead cells and a compensation matrix was given
  if (!is.null(ld_channel) & !is.null(compensation)) {
    # Set default parameters for live gate
    live_defaults <- list("channels" = c(ld_channel, "FSC-A"),
                          "position" = c(FALSE, NA),
                          "percentile" = c(0.70, NA),
                          "twin.factor" = c(0.1, NA))
    # Add to or override defaults with user specified parameters
    if (is.null(live_gate)) {
      live_gate <- live_defaults
    } else if (is.list(live_gate)) {
      live_gate <- utils::modifyList(live_defaults, live_gate)
    } else {
      live_gate <- prepareGateParams(live_gate)
    }
  }

  # Remove margin events and doublets
  ff_m <- PeacoQC::RemoveMargins(ff, c("FSC-A", "SSC-A", all_channels))
  ff_d <- PeacoQC::RemoveDoublets(ff_m, nmad = nmad)

  # Apply debris gate
  ff_g <- do.call(flowDensity::flowDensity, c(list("obj" = ff_d), debris_gate))
  ff_g <- flowDensity::getflowFrame(ff_g)

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
                             "logicle" = flowCore::estimateLogicle(ff_c, channels = ld_channel),
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
    ff_l <- do.call(flowDensity::flowDensity, c(list("obj" = ff_t), live_gate))
    ff_l <- flowDensity::getflowFrame(ff_l)
    # o = flowDensity::deGate(ff_t, "BUV496-A", tinypeak.removal = 1/20, spar = 0.9)
    # print(o)
    # rect = flowCore::rectangleGate(.gate = list("BUV496-A"= c(-0.5, o),
    #                                             "FSC-A" = c(0, 300000)))
    # filt = flowCore::filter(ff_t, rect)
    # ff_l = flowCore::Subset(ff_t, filt)
  }

  # Create plot for results of doublet removal
  p1 <- plotBeforeAfter(ff_m, ff_d, "FSC-A", "FSC-H", 5000)

  # Create plot for results of debris removal
  p2 <- plotBeforeAfter(ff_d, ff_g, "FSC-A", "SSC-A", 5000)

  # Create plot for results of live/dead cell removal
  p3 <- plotBeforeAfter(ff_t, ff_l, ld_channel, "FSC-A", 5000)

  # Arrange above plots
  gridExtra::grid.arrange(p1, p2, p3,
                          nrow = 2,
                          ncol = 2,
                          layout_matrix = rbind(c(1,1,2,2),
                                                c(1,1,2,2),
                                                c(NA,3,3,NA),
                                                c(NA,3,3,NA)),
                          top = file)
}


#' previewPreprocessing.GatingSet
#'
#' ...
#'
#' @inheritParams doPreprocessing.GatingSet
#' @param sample_ind A character string or integer specifying which sample to examine;
#' default is \code{1}
#'
#' @importFrom patchwork wrap_plots
#' @importFrom flowStats gate_singlet
#'
#' @return A plot
#' @export
previewPreprocessing.GatingSet <- function(input,
                                           sample_ind = 1,
                                           compensation = NULL,
                                           transformation = c("logicle", "arcsinh", "linear"),
                                           ld_channel = NULL,
                                           debris_args = list(),
                                           singlet_args = list(),
                                           live_args = list()) {
  # Get copy of first sample
  ff <- flowWorkspace::cytoframe_to_flowFrame(flowWorkspace::gs_pop_get_data(input)[[sample_ind]])

  # If a compensation matrix was given
  if (!is.null(compensation)) {
    # Apply compensation matrix
    ff <- flowCore::compensate(ff, comp_mat)
  }

  # Apply transformation if it was given/desired
  if (!is.null(transformation)) {
    all_channels <- names(flowWorkspace::markernames(input))
    if(is.character(transformation)) {
      transformation <- switch(transformation,
                               "logicle" = flowCore::estimateLogicle(ff, channels = all_channels),
                               "arcsinh" = {
                                 res <- flowCore::arcsinhTransform(b = 5)
                                 res <- flowCore::transformList(all_channels, res)
                                 res},
                               "linear" = {
                                 res <- flowCore::linearTransform()
                                 res <- flowCore::transformList(all_channels, res)
                                 res
                               })
    }
    ff <- flowCore::transform(ff, transformation)
  }

  # Add any additional options chosen by user, and overwrite defaults if needed
  debris_options <- utils::modifyList(list(channel = "FSC-A"), debris_args)
  live_options <- utils::modifyList(list(channel = ld_channel, positive = FALSE), live_args)

  # Generate gates
  g1 <- do.call(openCyto::gate_mindensity, c(list(fr = ff), debris_options))
  g2 <- do.call(flowStats::gate_singlet, c(list(x = ff), singlet_args))
  g3 <- do.call(openCyto::gate_mindensity, c(list(fr = ff), live_options))

  # Make gates
  # g1 <- openCyto::gate_mindensity(ff, channel = "FSC-A")
  # g2 <- flowStats::gate_singlet(ff)
  # g3 <- openCyto::gate_mindensity(ff, channel = ld_channel, positive = FALSE)

  # # Make each plot
  # p1 <- ggcyto::as.ggplot(ggcyto::ggcyto(ff, aes(x = `FSC-A`, y = `SSC-A`)) +
  #                           ggplot2::geom_hex(bins = 200))
  # p2 <- ggcyto::as.ggplot(ggcyto::ggcyto(ff, aes(x = `FSC-A`, y = `FSC-H`)) +
  #                           ggplot2::geom_hex(bins = 200))
  # p3 <- ggcyto::as.ggplot(ggcyto::ggcyto(ff, aes(x = !!enquo(ld_channel), y = `FSC-A`)) +
  #                           ggplot2::geom_hex(bins = 200))

  ranges <- Biobase::pData(flowCore::parameters(ff)) %>%
    dplyr::filter(name == "FSC-A")
  lwr <- ranges$minRange
  upr <- ranges$maxRange
  # should filter out instead

  p1 <- ggcyto::as.ggplot(ggcyto::autoplot(ff, "FSC-A", "SSC-A", bins = 200) +
                            ggplot2::xlim(c(lwr, upr)) +
                            ggcyto::geom_gate(g1))
  ff <- flowCore::Subset(ff, g1)

  p2 <- ggcyto::as.ggplot(ggcyto::autoplot(ff, "FSC-A", "FSC-H", bins = 200) +
                            ggplot2::xlim(c(lwr, upr)) +
                            ggcyto::geom_gate(g2))
  ff <- flowCore::Subset(ff, g2)

  p3 <- ggcyto::as.ggplot(ggcyto::autoplot(ff, ld_channel, "FSC-A", bins = 200) +
                            ggcyto::geom_gate(g3))

  # Combine plots into one image
  patchwork::wrap_plots(p1, p2, p3, ncol = 1)

}
