#' generateGatingTable
#'
#' @param gs A `GatingSet`
#' @param collapse_data A `boolean`, indicating if samples should be collapsed
#' into one
#' @param ld_stain If the experiment contains a live/dead stain, a `character`
#' indicating which channel
#'
#' @return A data.table, to be used for creation of a `gatingTemplate`
#' @seealso [openCyto::gatingTemplate()]
#'
#' @export
generateGatingTable <- function(gs, collapse_data = FALSE, ld_stain = NULL) {
  # Get number of samples in experiment
  num_samples <- length(gs)

  group_by <- NA
  if(collapse_data) {
    group_by <- num_samples
  }

  # Make default gating table
  gt_table <- data.table(alias = c("nonMargins", "nonDebris", "singlets"),
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
    live_gate <- data.table(alias = "live",
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
#' @return A \code{ggplot} drawn with \code{ggcyto}
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

#' editGateManual
#'
#' A function enabling the user to redraw gates for all or a subset of a
#' \code{GatingSet}. !!! not working with current versions of ggcyto and ggplot2
#'
#' @param gs A \code{GatingSet} whose gates will be manually edited
#' @param node A \code{character} string specifying which node to redraw
#' @param dims A \code{character} vector indicating which channels to draw the new
#' gate on. By default, it is set to the channels used in the original gate. If it
#' is of length-1, a histogram will be displayed for the user to draw on. If
#' length-2, a scatterplot will be displayed instead.
#' @param ref_sample A \code{numeric} or \code{character}, specifying the index
#' or name of the sample within the \code{GatingSet} that you would like to
#' draw the new gate on. NOTE: if \code{sample_ids} is specified by the user,
#' \code{ref_sample} should be found within the vector \code{sample_ids}.
#'  (e.g. if I have a GatingSet
#'  of length 10, and set \code{sample_ids = 6:10}, then \code{ref_sample} should
#'  be a number between 6 and 10.)
#' @param sample_ids Optional; if you would like to redraw this gate for only a subset
#' of your samples, a vector of sample indices or names.
#'
#' @return A \code{GatingSet} whose gates have been redrawn according to user input.
#' @export
editGateManual <- function(gs, node, dims = NULL, ref_sample = 1, sample_ids = NULL) {

  # Subset GatingSet if sample_ids given
  if (!is.null(sample_ids)) {
    sn <- flowWorkspace::sampleNames(gs)
    gs <- gs[sample_ids]

    viable_ids <- match(flowWorkspace::sampleNames(gs), sn)

    if (length(viable_ids == 1)) {
      ref_sample <- 1
    }
    else if (!any(ref_sample %in% viable_ids)) {
      warning(paste0("\n `ref_sample=", ref_sample,
                     "` is not found in the given `sample_ids`, should be one of: ",
                     paste(viable_ids, collapse=","), ". \n \n Using `ref_sample=",
                     viable_ids[1], "` by default."),
              immediate. = TRUE)

      ref_sample <- 1
    }
    else {
      ref_sample <- match(sn, flowWorkspace::sampleNames(gs))[ref_sample]
    }
  }

  # Get default dimensions if none were given
  if (is.null(dims)) {
    dims <- flowWorkspace::gs_pop_get_gate(gs, node)[[1]]
    dims <- names(methods::slot(dims, "parameters"))
  }


  # Get immediate parent of node to be redrawn
  node_parent <- flowWorkspace::gs_pop_get_parent(gs, node, path = "full")


  # Open R Shiny window to redraw gate
  flowGate::gs_gate_interactive(gs,
                                filterId = "newGate",
                                sample = ref_sample, # sample to draw the gate on
                                subset = node_parent,
                                dims = dims,
                                regate = FALSE, # if true, redrawn gate and children are deleted
                                overlayGates = node # overlay gate that is being redrawn
  )

  # Get new gate that was just drawn
  gate <- flowWorkspace::gs_pop_get_gate(gs, "newGate")

  # Set existing gate to the new one
  flowWorkspace::gs_pop_set_gate(gs, node, gate)

  # Recompute descendants of the redrawn gate
  flowWorkspace::recompute(gs, node)

  # Delete temporary gate
  flowWorkspace::gs_pop_remove(gs, "newGate")
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
