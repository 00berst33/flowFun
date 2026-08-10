#' generateGatingTable
#'
#' Generates a default gating template for preprocessing
#'
#' @param gs A `GatingSet`
#' @param collapse_data A `boolean`, indicating if samples should be collapsed
#' into one when drawing gates. Default is `TRUE`
#' @param ld_stain If the experiment contains a live/dead stain, a `character`
#' indicating which channel
#'
#' @details
#' The package `openCyto` provides functions to automatically draw gates for cytometry
#' data. To work, it requires a \code{\link[openCyto:gatingTemplate]{gatingTemplate}}
#' object specifying basic information about the gates to be drawn.
#'
#' This function generates a simple gating scheme for use with `openCyto`
#' containing the following steps:
#'  - Remove margin events
#'  - Remove debris
#'  - Remove doublets
#'  - Remove dead cells (if a viability stain was used in the panel)
#'
#' The argument `collapse_data` affects how gates are drawn. If `TRUE` (default),
#' all samples are aggregated into one, and gates are drawn using said aggregated sample.
#' If `FALSE`, gates are instead drawn separately for each sample.
#'
#' @return A `data.table`, to be used for creation of a `gatingTemplate`
#' @seealso \code{\link[openCyto:gatingTemplate]{gatingTemplate()}},
#' \code{\link[openCyto:gatingTemplate]{gt_gating()}}
#'
#' @export
generateGatingTable <- function(gs, collapse_data = TRUE, ld_stain = NULL) {
  # Get number of samples in experiment
  num_samples <- length(gs)

  # Set values to collapse samples, if desired
  group_by <- NA
  if(collapse_data) {
    group_by <- num_samples
  }

  # Make default gating table
  gt_table <- data.table::data.table(alias = c("nonMargins", "nonDebris", "singlets"),
                                     pop = c("+", "+", "+"),
                                     parent = c("root", "nonMargins", "nonDebris"),
                                     dims = c("FSC-A,SSC-A", "FSC-A", "FSC-A,FSC-H"),
                                     gating_method = c("boundary", "gate_mindensity", "singletGate"),
                                     gating_args = c("min=c(0,0),max=c(262143,262143)", NA, NA),
                                     collapseDataForGating = rep(collapse_data, 3),
                                     groupBy = as.integer(rep(group_by, 3)),
                                     preprocessing_method = c(NA, NA, NA),
                                     preprocessing_args = c(NA, NA, NA))

  # If a L/D stain was given
  if(!is.null(ld_stain)) {
    # Make gate for live cells
    live_gate <- data.table::data.table(alias = "live",
                                        pop = "-",
                                        parent = "singlets",
                                        dims = "BUV496-A",
                                        gating_method = "gate_mindensity",
                                        gating_args = NA,
                                        collapseDataForGating = collapse_data,
                                        groupBy = as.integer(group_by),
                                        preprocessing_method = NA,
                                        preprocessing_args = NA)

    # Add to gating table
    gt_table <- rbind(gt_table, live_gate)
  }

  return(gt_table)
}


#' plotAllSamples
#'
#' Plot chosen gate for all samples
#'
#' @param gs A \code{GatingSet}
#' @param xdim The channel to plot on the x-axis
#' @param ydim The channel to plot on the y-axis
#' @param subset The population whose cells should be displayed on the plot
#' @param node The gate to be visualized
#'
#' @details
#' Intended for checking the results of newly added gates, typically after a
#' call to \code{\link[openCyto:gt_gating]{gt_gating()}}. Note that to use this
#' function, `gs` must have at least one gate.
#'
#' @return A \code{ggplot} drawn with \code{ggcyto}
#'
#' @examples
#' # load in data
#' files <- system.file("extdata", "samples", package = "flowFunData")
#' cs <- flowWorkspace::load_cytoset_from_fcs(path=dir, pattern=".fcs", truncate_max_range = TRUE)
#'
#' # Make GatingSet
#' gs <- flowWorkspace::GatingSet(cs)
#'
#' # Generate default gating scheme
#' gt_table <- generateGatingTable(gs, collapse_data = TRUE)
#' gt <- openCyto::gatingTemplate(gt_table)
#'
#' # Apply gating scheme
#'   # not run: assume transformations and compensations were performed as needed
#' openCyto::gt_gating(gt, gs)
#'
#' ## Plotting results
#' # View nonDebris gate
#' plotAllSamples(gs, "FSC-A", "SSC-A", "nonMargins", "nonDebris")
#'
#' # View singlets gate
#' plotAllSamples(gs, "FSC-A", "FSC-H", "nonDebris", "singlets")
#'
#' @export
plotAllSamples <- function(gs, xdim, ydim, subset, node) {
  xdim <- rlang::enquo(xdim)
  ydim <- rlang::enquo(ydim)
  p <- ggcyto::ggcyto(gs, mapping = ggplot2::aes(x = !!xdim, y = !!ydim), subset = subset) +
    ggplot2::geom_hex(bins = 50) +
    ggcyto::geom_gate(node) +
    ggplot2::theme(text = ggplot2::element_text(size = 10)) +
    ggcyto::geom_stats(size = 4)

  return(p)
}


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
