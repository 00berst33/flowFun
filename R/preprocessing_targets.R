tar_option_set(
  packages = c("flowCore", "PeacoQC", "Biobase", "ggplot2"), # Packages that your targets need for their tasks.
  format = "qs"
  # Set other options as needed.
)

#' @keywords internal
#' @export
getPanel <- function(ff) {
  ff_df <- Biobase::pData(flowCore::parameters(ff))
  used_ind <- which(!is.na(ff_df$desc))
  panel_vec <- ff_df$name[used_ind]
  names(panel_vec) <- ff_df$desc[used_ind]
  return(panel_vec)
}

#' @keywords internal
#' @export
getTransform <- function(ff, transformation, transformation_type, panel) {
  panel <- as.vector(panel)

  if (is.null(transformation)) {
    transformation <- switch(transformation_type,
                             "logicle" = flowCore::estimateLogicle(ff, channels = panel),
                             "arcsinh" = {
                               transformation <- flowCore::arcsinhTransform(b = 5)
                               transformation <- flowCore::transformList(panel, tfun = transformation)
                             },
                             "none" = {
                               transformation <- flowCore::linearTransform()
                               transformation <- flowCore::transformList(panel, tfun = transformation)
                             },
                             "other" = stop("Must provide a value for `transformation` when `transformation_type` is 'other'."))
  }

  return(transformation)
}


drawGrobs <- function(p1, p2, p3) {
  # p <- list()
  # i <- 1
  #
  # # Create plot for results of doublet removal
  # p1 <- ggplot2::ggplotGrob(p1)
  # p[[i]] <- p1
  #
  # # Create plot for results of debris removal
  # p2 <- ggplot2::ggplotGrob(p2)
  # p[[i]] <- p2
  #
  # # Create plot for results of live/dead cell removal
  # if (!is.null(ld_channel)) {
  #   p3 <- ggplot2::ggplotGrob(p3)
  #   p[[i]] <- p3
  # }
  p <- as.list(p1, p2, p3)
  p <- lapply(p, function(i) {if (!is.null(p[[i]])) {p[[i]] <- ggplot2::ggplotGrob(p[[i]])}})

  # Arrange above plots
  p <- gridExtra::marrangeGrob(p,
                          nrow = 2,
                          ncol = 2,
                          layout_matrix = rbind(c(1,1,2,2),
                                                c(1,1,2,2),
                                                c(NA,3,3,NA),
                                                c(NA,3,3,NA)),
                          top = file)
  return(p)
}

#' @export
flowPreprocessing <- function(file, compensation = NULL, transformation = NULL,
                              transformation_type = c("logicle", "arcsinh", "other", "none"),
                              debris_gate = NULL, ld_channel = NULL, live_gate = NULL, nmad = 4) {

  list(
    targets::tar_target_raw("read_file",
                            substitute(flowCore::read.FCS(file, truncate_max_range = FALSE)),
                            format = "qs"),
    targets::tar_target_raw("panel", substitute(getPanel(read_file)), format = "qs"),
    targets::tar_target_raw("margins_removed",
                            substitute(PeacoQC::RemoveMargins(read_file, c("FSC-A", panel))),
                            format = "qs"),
    targets::tar_target_raw("doublets_removed",
                            substitute(PeacoQC::RemoveDoublets(margins_removed, nmad = nmad)),
                            format = "qs"),
    targets::tar_target_raw("rem_doublets_plot",
                            substitute(plotBeforeAfter(margins_removed, doublets_removed, "FSC-A", "FSC-H", 5000)),
                            format = "qs"),
    targets::tar_target_raw("debris_removed",
                            substitute(flowDensity::getflowFrame(flowDensity::flowDensity(
                              doublets_removed,
                              channels = c("FSC-A", "SSC-A"),
                              position = c(TRUE, NA),
                              percentile = c(0.85, NA)))),
                            cue = tarchetypes::tar_cue_skip(is.null(debris_gate))),
    targets::tar_target_raw("rem_debris_plot",
                            substitute(plotBeforeAfter(doublets_removed, debris_removed, "FSC-A", "SSC-A", 5000))),
    targets::tar_target_raw("compensation_applied",
                            substitute(flowCore::compensate(debris_removed,
                                                            spillover = compensation)),
                            cue = tarchetypes::tar_cue_skip(is.null(compensation))),
    targets::tar_target_raw("transform_applied",
                            substitute(flowCore::transform(compensation_applied,
                                                           getTransform(compensation_applied,
                                                                        transformation,
                                                                        transformation_type,
                                                                        panel)))),
    targets::tar_target_raw("dead_removed",
                            substitute(flowDensity::getflowFrame(flowDensity::flowDensity(
                              transform_applied,
                              channels = c(ld_channel, "FSC-A"),
                              position = c(FALSE, NA),
                              percentile = c(0.70, NA),
                              twin.factor = c(0.1, NA)))),
                              cue = tarchetypes::tar_cue_skip(is.null(ld_channel) | is.null(compensation))),
    targets::tar_target_raw("rem_dead_plot",
                            substitute(plotBeforeAfter(transform_applied, dead_removed, ld_channel, "FSC-A", 5000)),
                            cue = tarchetypes::tar_cue_skip(is.null(ld_channel) | is.null(compensation))),
    targets::tar_target_raw("qc_applied",
                            substitute(flowCut::flowCut(f = dead_removed,
                                                        FileID = make.names(file),
                                                        Plot = "None")[[1]]))
    # targets::tar_target_raw("all_rem_plots",
    #                         substitute(drawGrobs(rem_doublets_plot, rem_debris_plot, rem_dead_plot)))
  )
}

# run_pi <- function() {
#   library(targets)
#   # library(flowFun)
#
#   # Set target options:
#   targets::tar_option_set(
#     packages = c("flowCore", "PeacoQC", "Biobase", "ggplot2"), # Packages that your targets need for their tasks.
#     format = "qs"
#     # Set other options as needed.
#   )
#   s
#   # Replace the target list below with your own:
#   file <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw/Ab_PHA_Ctrl AWB4.fcs"
#   # comp <- system.file("extdata", "compensation_matrix.csv", package = "flowFun")
#
#   flowPreprocessing(file, compensation = comp, transformation_type = "logicle")
# }
