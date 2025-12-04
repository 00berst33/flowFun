#' getMetaclusterMFIs
#'
#' @param input A data frame or table, with column \code{"Metacluster"}.
#' @param cols_to_use A vector specifying which columns to calculate medians for.
#' Default is all columns used for clustering.
#'
#' @return A data frame, where rows are metaclusters and columns are channels.
#'
#' @export
getMetaclusterMFIs <- function(input, cols_to_use = NULL) {
  Metacluster <- NULL

  if (is.null(cols_to_use)) {
    cols_to_use <- attr(input, "clustered")
  }

  medians <- input %>%
    tidytable::group_by(Metacluster) %>%
    tidytable::summarise(tidytable::across(.cols = cols_to_use,
                                           .fns = stats::median,
                                           .drop = "keep")) %>%
    tidytable::select(c(Metacluster, cols_to_use)) %>%
    data.frame(row.names = "Metacluster",
               check.names = FALSE)

  rownames(medians) <- basename(rownames(medians))

  return(medians)
}

#' getClusterMFIs
#'
#' @param input A data frame or table, with column \code{"Cluster"}.
#' @param cols_to_use A vector specifying which columns to calculate medians for.
#' Default is all columns used for clustering.
#'
#' @return A data frame, where rows are clusters and columns are channels.
#'
#' @export
getClusterMFIs <- function(input, cols_to_use = NULL) {
  Cluster <- NULL

  if (is.null(cols_to_use)) {
    cols_to_use <- attr(input, "clustered")
  }

  medians <- input %>%
    tidytable::group_by(Cluster) %>%
    tidytable::summarise(tidytable::across(.cols = cols_to_use,
                                           .fns = stats::median,
                                           .drop = "keep")) %>%
    tidytable::select(c(Cluster, cols_to_use)) %>%
    data.frame(row.names = "Cluster",
               check.names = FALSE)

  rownames(medians) <- basename(rownames(medians))

  return(medians)
}

#' plotMetaclusterMFIs
#'
#' Plot a heatmap of metacluster MFIs.
#'
#' @param input A FlowSOM object or data table.
#' @param cols_to_use A character vector specifying which markers you
#' would like to include in the plot. Default is the markers used for clustering.
#' @param ... Additional parameters to pass to [ComplexHeatmap::Heatmap()].
#'
#' @return A heatmap of MFIs made with \code{ComplexHeatmap::Heatmap()},
#' where each column is a marker of interest, and each row is a metacluster.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
#' fsom <- FlowSOM::FlowSOM(file,
#'                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#'                          nClus = 10,
#'                          seed = 42,
#'                          xdim = 6,
#'                          ydim = 6)
#'
#' plotMetaclusterMFIs(fsom)
plotMetaclusterMFIs <- function(input, cols_to_use = NULL, ...) {
  result <- UseMethod("plotMetaclusterMFIs")
  return(result)
}

#' plotMetaclusterMFIs.FlowSOM
#'
#' @keywords internal
#' @export
plotMetaclusterMFIs.FlowSOM = function(input, cols_to_use = NULL, ...) {
  if (is.null(cols_to_use)) {
    cols_to_use <- input$map$colsUsed
  }

  # Get metacluster MFIs for each marker/channel of interest
  mfi_mat <- FlowSOM::GetMetaclusterMFIs(fsom = input, colsUsed = FALSE, prettyColnames = FALSE)
  mfi_mat <- mfi_mat[, which(colnames(mfi_mat) %in% FlowSOM::GetChannels(input, cols_to_use))]
  colnames(mfi_mat) <- input$prettyColnames[colnames(mfi_mat)]
  rownames(mfi_mat) <- levels(FlowSOM::GetMetaclusters(input))

  # Set default heatmap options
  default_options <- list(border = TRUE,
                          show_row_names = TRUE,
                          show_column_dend = FALSE,
                          row_names_gp = grid::gpar(fontsize = 7),
                          column_names_gp = grid::gpar(fontsize = 7),
                          heatmap_legend_param = list(
                           title = "Expression",
                           title_gp = grid::gpar(fontsize = 9),
                           labels_gp = grid::gpar(fontsize = 8)))

  # Add any additional options chosen by user, and overwrite defaults if needed
  additional_options <- list(...)
  heatmap_options <- utils::modifyList(default_options, additional_options)

  # Generate heatmap
  mfi_heatmap <- do.call(ComplexHeatmap::Heatmap, c(list(matrix = as.matrix(mfi_mat)), heatmap_options))

  return(mfi_heatmap)
}

#' plotMetaclusterMFIs.data.frame
#'
#' @keywords internal
#' @export
plotMetaclusterMFIs.data.frame <- function(input, cols_to_use = NULL, ...) {
  Metacluster <- NULL

  if (methods::is(input, "data.table") & is.null(cols_to_use)) {
    cols_to_use <- attributes(input)$clustered
  }

  # Get metacluster MFIs for each marker/channel of interest
  mfi_mat <- input %>%
    tidytable::summarise(tidytable::across(.cols = cols_to_use, # column
                                           .fns = stats::median,
                                           .drop = "keep"),
                         .by = Metacluster)
  mfi_mat <- as.matrix(mfi_mat, rownames = "Metacluster")

  colnames(mfi_mat) <- sapply(colnames(mfi_mat), getPrettyColNames, input = input)

  # Set default heatmap options
  default_options <- list(border = TRUE,
                          show_row_names = TRUE,
                          show_column_dend = FALSE,
                          row_names_gp = grid::gpar(fontsize = 7),
                          column_names_gp = grid::gpar(fontsize = 7),
                          heatmap_legend_param = list(
                            title = "Expression",
                            title_gp = grid::gpar(fontsize = 9),
                            labels_gp = grid::gpar(fontsize = 8)))

  # Add any additional options chosen by user, and overwrite defaults if needed
  additional_options <- list(...)
  heatmap_options <- utils::modifyList(default_options, additional_options)

  # Generate heatmap
  mfi_heatmap <- do.call(ComplexHeatmap::Heatmap, c(list(matrix = mfi_mat), heatmap_options))

  return(mfi_heatmap)
}

getPrettyColNames <- function(col, input) {
  temp_col <- which(colnames(input) == col)
  colnames(input)[temp_col] <- paste0(col, " <", attr(input[[col]], "marker"), ">")
}

#' plotClusterMFIs
#'
#' Plot a heatmap of cluster MFIs.
#'
#' @param input The FlowSOM object or data table to plot.
#' @param cols_to_use A character vector specifying which markers you
#' would like to include in the plot. Default is the markers
#' used for clustering.
#' @param metaclusters A vector specifying which metaclusters whose clusters
#' you would like to include in the plot. If NULL (default), all clusters are included.
#' @param ... Additional parameters to pass to [ComplexHeatmap::Heatmap()].
#'
#' @return A heatmap of MFIs made with \code{ComplexHeatmap::Heatmap()}, where each
#' column is a marker of interest, and each row is a cluster.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
#' fsom <- FlowSOM::FlowSOM(file,
#'                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#'                          nClus = 10,
#'                          seed = 42,
#'                          xdim = 6,
#'                          ydim = 6)
#'
#' # Plot all clusters
#' plotClusterMFIs(fsom)
#'
#' # Plot only clusters belonging to metaclusters 9 and 10
#' plotClusterMFIs(fsom, metaclusters = c(9, 10))
plotClusterMFIs <- function(input, cols_to_use = NULL, metaclusters = NULL, ...) {
  results <- UseMethod("plotClusterMFIs")
  return(results)
}

#' @keywords internal
#' @export
plotClusterMFIs.FlowSOM <- function(input, cols_to_use = input$map$colsUsed,
                                    metaclusters = NULL, ...) {

  # Get cluster MFIs for each marker/channel of interest
  mfi_mat <- FlowSOM::GetClusterMFIs(input, colsUsed = FALSE, prettyColnames = FALSE)
  mfi_mat <- mfi_mat[, which(colnames(mfi_mat) %in% FlowSOM::GetChannels(input, cols_to_use))]

  # Plot only clusters belonging to given metaclusters
  if (!is.null(metaclusters)) {
    ind <- input$metaclustering %in% metaclusters
    if (!(TRUE %in% ind)) {
      stop("No clusters found in the given metacluster(s). Check that either your indices are within bounds, or names have no typos.")
    }
    mfi_mat <- mfi_mat[ind, ]
  }

  colnames(mfi_mat) <- input$prettyColnames[colnames(mfi_mat)]

  # Set default heatmap options
  default_options <- list(border = TRUE,
                         show_row_names = TRUE,
                         show_column_dend = FALSE,
                         row_names_gp = grid::gpar(fontsize = 4),
                         column_names_gp = grid::gpar(fontsize = 6),
                         width = grid::unit(0.5, "npc"),
                         heatmap_legend_param = list(
                           title = "Expression",
                           title_gp = grid::gpar(fontsize = 9),
                           labels_gp = grid::gpar(fontsize = 8)))

  # Add any additional options chosen by user, overwrite defaults if necessary
  additional_options <- list(...)
  heatmap_options <- utils::modifyList(default_options, additional_options)

  # Generate heatmap
  mfi_heatmap <- do.call(ComplexHeatmap::Heatmap, c(list(matrix = as.matrix(mfi_mat)), heatmap_options))

  return(mfi_heatmap)
}

#' @keywords internal
#' @export
plotClusterMFIs.data.frame <- function(input, cols_to_use = NULL,
                                       metaclusters = NULL, ...) {
  Metacluster <- Cluster <- NULL

  if (methods::is(input, "data.table") & is.null(cols_to_use)) {
    cols_to_use <- attributes(input)$clustered
  }

  # Get cluster MFIs for each marker/channel of interest
  mfi_mat <- input %>%
    tidytable::summarise(tidytable::across(.cols = cols_to_use,
                                           .fns = stats::median,
                                           .drop = "keep"),
                         .by = Cluster)

  # If metacluster(s) of interest were specified
  if (!is.null(metaclusters)) {
    # Get only clusters belonging to given metaclusters
    clust_to_keep <- input %>%
      tidytable::filter(Metacluster %in% metaclusters) %>%
      tidytable::pull(Cluster) %>%
      unique()

    # If no clusters belong to the given metacluster(s)
    if (length(clust_to_keep) < 1) {
      stop("No clusters found in the given metacluster(s). Check that either your indices are within bounds, or names have no typos.")
    }

    # Subset MFI matrix to relevant clusters
    mfi_mat <- mfi_mat %>%
      tidytable::filter(Cluster %in% clust_to_keep)
  }

  # Turn results into a matrix and rename columns
  mfi_mat <- as.matrix(mfi_mat, rownames = "Cluster")
  colnames(mfi_mat) <- sapply(colnames(mfi_mat), getPrettyColNames, input = input)

  # Set default heatmap options
  default_options <- list(border = TRUE,
                          show_row_names = TRUE,
                          show_column_dend = FALSE,
                          row_names_gp = grid::gpar(fontsize = 4),
                          column_names_gp = grid::gpar(fontsize = 6),
                          width = grid::unit(0.5, "npc"),
                          heatmap_legend_param = list(
                            title = "Expression",
                            title_gp = grid::gpar(fontsize = 9),
                            labels_gp = grid::gpar(fontsize = 8)))

  # Add any additional options chosen by user, overwrite defaults if necessary
  additional_options <- list(...)
  heatmap_options <- utils::modifyList(default_options, additional_options)

  # Generate heatmap
  mfi_heatmap <- do.call(ComplexHeatmap::Heatmap, c(list(matrix = as.matrix(mfi_mat)), heatmap_options))

  return(mfi_heatmap)
}

#' annotateMFIHeatmap
#'
#' Annotate an unmerged heatmap of MFIs with each metacluster's final assignment.
#'
#' @param merged_input The final table or FlowSOM object, whose metaclusters have
#' been merged and named via either [editTableMetaclusters()] or
#' [FlowSOM::UpdateMetaclusters()].
#' @param original_input If \code{merged_input} is a FlowSOM object, the initial
#' FlowSOM object before merging should be provided here.
#' @param cols_to_use A character vector specifying which markers you
#' would like to include in the plot. Default is the markers used for clustering.
#' @param ... Additional parameters to pass to [ComplexHeatmap::Heatmap()].
#'
#' @details
#' Additional details...
#'
#' @return An annotated heatmap generated by [ComplexHeatmap::Heatmap()].
#'
#' @export
annotateMFIHeatmap <- function(merged_input, original_input = NULL,
                               cols_to_use = NULL, ...) {
  result <- UseMethod("annotateMFIHeatmap")
  return(result)
}

#' @keywords internal
#' @export
annotateMFIHeatmap.FlowSOM <- function(merged_input, original_input, cols_to_use = NULL, ...) {
  if (methods::is(merged_input, "data.table") & is.null(cols_to_use)) {
    cols_to_use <- attributes(merged_input)$clustered
  }

  # Generate heatmap
  heatmap <- plotMetaclusterMFIs(merged_input, cols_to_use)

  # Get final metacluster assignments
  metaclusters <- levels(original_input$metaclustering)
  new_assignments <- c()
  for (metacluster in metaclusters) {
    cell_ind <- which(original_input$metaclustering == metacluster)[1]

    new_meta <- as.vector(merged_input$metaclustering[cell_ind])

    new_assignments <- c(new_assignments, new_meta)
  }

  # Define colors for annotation
  row_colors <- viridis::viridis(length(metaclusters), option = "H")[as.numeric(factor(new_assignments))]
  names(row_colors) <- new_assignments

  # Add annotation to heatmap
  anno <- ComplexHeatmap::rowAnnotation(Metacluster = new_assignments,
                                        show_annotation_name = FALSE,
                                        col = list(Metacluster = row_colors))
  return(heatmap + anno)
}

#' @keywords internal
#' @export
annotateMFIHeatmap.data.frame <- function(merged_input, cols_to_use = NULL, ...) {
  Meta_original <- NULL

  # Get columns to use if necessary
  if (methods::is(merged_input, "data.table") & is.null(cols_to_use)) {
    cols_to_use <- attributes(merged_input)$clustered
  }

  # Get metacluster MFIs for each marker/channel of interest
  mfi_mat <- merged_input %>%
    tidytable::summarise(tidytable::across(.cols = cols_to_use, # column
                                           .fns = stats::median,
                                           .drop = "keep"),
                         .by = Meta_original)
  mfi_mat <- as.matrix(mfi_mat, rownames = "Meta_original")
  colnames(mfi_mat) <- sapply(colnames(mfi_mat), getPrettyColNames, input = merged_input)

  # Set default heatmap options
  default_options <- list(border = TRUE,
                          show_row_names = TRUE,
                          show_column_dend = FALSE,
                          row_names_gp = grid::gpar(fontsize = 7),
                          column_names_gp = grid::gpar(fontsize = 7),
                          heatmap_legend_param = list(
                            title = "Expression",
                            title_gp = grid::gpar(fontsize = 9),
                            labels_gp = grid::gpar(fontsize = 8)))

  # Add any additional options chosen by user, and overwrite defaults if needed
  additional_options <- list(...)
  heatmap_options <- utils::modifyList(default_options, additional_options)

  # Generate heatmap
  heatmap <- do.call(ComplexHeatmap::Heatmap, c(list(matrix = mfi_mat), heatmap_options))

  # Get final metacluster assignments
  metaclusters <- levels(merged_input$Meta_original)
  new_assignments <- c()
  for (metacluster in metaclusters) {
    cell_ind <- which(merged_input$Meta_original == metacluster)[1]

    new_meta <- as.vector(merged_input$Metacluster[cell_ind])

    new_assignments <- c(new_assignments, new_meta)
  }

  # Define colors for annotation
  row_colors <- viridis::viridis(length(metaclusters), option = "H")[as.numeric(factor(new_assignments))]
  names(row_colors) <- new_assignments

  # Add annotation to heatmap
  anno <- ComplexHeatmap::rowAnnotation(Metacluster = new_assignments,
                                        show_annotation_name = FALSE,
                                        col = list(Metacluster = row_colors))
  return(heatmap + anno)
}

#' searchByExpression
#'
#' Search MST for nodes with certain marker expressions.
#'
#' @param fsom A FlowSOM object as generated by the FlowSOM function.
#' @param levels A list specifying channels of interest and whether they are
#' high or low. Note that searching for intermediate expression is not supported.
#' See examples for the proper format of this variable.
#'
#' @return An MST plot, where nodes with the given expression levels are
#' highlighted and labeled as "True".
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
#' fsom <- FlowSOM::FlowSOM(file,
#'                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#'                          nClus = 10,
#'                          seed = 42,
#'                          xdim = 6,
#'                          ydim = 6)
#'
#' searchByExpression(fsom,
#'                    levels = c("BUV661-A" = "high",
#'                               "APC-Cy7-A" = "low"))
#'
#' searchByExpression(fsom,
#'                    levels = c("BUV805-A" = "low",
#'                               "APC-Cy7-A" = "high",
#'                               "Alexa Fluor 700-A" = "high"))
searchByExpression = function(fsom, levels) {
  query <- levels
  query_res <- FlowSOM::QueryStarPlot(fsom, query, equalNodeSize = TRUE, plot = FALSE)
  cell_types <- factor(rep("False", fsom$map$nNodes),
                      levels = c("False", "True"))
  cell_types[query_res$selected] <- "True"
  FlowSOM::PlotStars(fsom, backgroundValues <- cell_types)
}

#' plotLabeled2DScatter
#'
#' Plot a 2D scatterplot colored by metacluster.
#'
#' @param input A FlowSOM object or data table.
#' @param channelpair A vector of two channel names.
#' @param clusters A vector specifying clusters of interest.
#' @param metaclusters A vector specifying metaclusters of interest.
#' @param labels Logical, should labels be plotted for each cluster center. Default
#' is \code{TRUE}.
#' @param colors Optional, a vector of colors for each metacluster.
#'
#' @details
#' It is only necessary to define one of the arguments \code{clusters} and
#' \code{metaclusters}, but the user may use both if they wish to plot a particular
#' subset of their data.
#'
#' This function draws a 2D scatterplot with the given channels on each axis.
#' Cells are colored by metacluster, and if \code{labels = TRUE}, a label is
#' plotted at each cluster center. The number contained in these labels indicates
#' the cluster's number, and the color of the label indicates which metacluster it
#' belongs to.
#'
#' @return A plot drawn with ggplot2.
#'
#' @export
#'
#' @examples
#' # Read in data
#' file <- system.file("extdata", "fsom_table_init.rds", package = "flowFun")
#' fsom_dt <- readRDS(file)
#'
#' # Generate plot
#' plotLabeled2DScatter(fsom_dt,
#'                      channelpair = c("APC-Cy7-A", "BUV615-P-A"),
#'                      metaclusters = c(17, 19))
#'
#' # Generate same plot with no labels
#' plotLabeled2DScatter(fsom_dt,
#'                      channelpair = c("APC-Cy7-A", "BUV615-P-A"),
#'                      metaclusters = c(17, 19),
#'                      labels = FALSE)
#'
#' # Specify both metaclusters and clusters
#' plotLabeled2DScatter(fsom_dt,
#'                      channelpair = c("APC-Cy7-A", "BUV563-A"),
#'                      clusters = c(66, 79, 91),
#'                      metaclusters = 16)
#'
#' # Specifying only clusters still colors by metacluster
#' plotLabeled2DScatter(fsom_dt,
#'                      channelpair = c("APC-Cy7-A", "BUV563-A"),
#'                      clusters = c(66, 79, 91))
plotLabeled2DScatter <- function(input, channelpair, clusters = NULL,
                                 metaclusters = NULL, labels = TRUE,
                                 colors = NULL) {
  result <- UseMethod("plotLabeled2DScatter")
  return(result)
}

### needs edits
# Note that the output of this function changes slightly depending on the combination
# of clusters and metaclusters the user is plotting. If the user plots metaclusters,
# with or without additional clusters defined in the \code{clusters} argument,
# then the output will be a single ggplot object, where each metacluster has a
# unique color. If the user plots only clusters, such that the \code{metaclusters}
# argument remains \code{NULL}, then a list of ggplot objects will be returned.
# Each plot will correspond to a cluster, and color only the cells belonging to
# said cluster. See examples for further explanation.
#
# The (meta)cluster number labels are included on every plot.
#
# file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
# fsom <- FlowSOM::FlowSOM(file,
#                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#                          nClus = 10,
#                          seed = 42,
#                          xdim = 6,
#                          ydim = 6)
#
# # Returns a single plot colored by metacluster
# plotLabeled2DScatter(fsom,
#                      channelpair = c("APC-Cy7-A", "BUV615-P-A"),
#                      metaclusters = c(1, 2, 7),
#                      label_type = "metacluster")
#
# # Also returns a single plot colored by metacluster
# plotLabeled2DScatter(fsom,
#                      channelpair = c("APC-Cy7-A", "BUV805-A"),
#                      clusters = c(25, 33),
#                      metaclusters = c(1, 4),
#                      label_type = "cluster")
#
# # Returns a list of two plots, once for each cluster
# p = plotLabeled2DScatter(fsom,
#                          channelpair = c("APC-Cy7-A", "BUV805-A"),
#                          clusters = c(4, 19),
#                          label_type = "cluster")
# plotLabeled2DScatter.FlowSOM <- function(input, channelpair, clusters = NULL,
#                                          metaclusters = NULL,
#                                          label_type = c("cluster", "metacluster")) {
#
#   if (is.null(clusters) && is.null(metaclusters)) {
#     stop("clusters and/or metaclusters must be a vector of indices.")
#
#   } else if (is.null(clusters) && !is.null(metaclusters)) {
#     # color by metacluster, plot colored labels with cluster numbers
#     all_clust <- which(input$metaclustering %in% metaclusters)
#
#   } else if (!is.null(clusters) && is.null(metaclusters)) {
#     # color by cluster
#     all_clust <- clusters
#
#   } else if (!is.null(clusters) && !is.null(metaclusters)) {
#     # color by metacluster, plot colored labels with cluster numbers
#     meta_nodes <- which(input$metaclustering %in% metaclusters)
#     all_clust <- unique(c(meta_nodes, clusters))
#
#   }
#
#   # edit !!!
#   # values in `metaclusters` and levels(input$metaclustering) can conflict
#   if (is.character(metaclusters)) {
#     metaclusters <- which(levels(input$metaclustering) %in% metaclusters)
#   } else if (is.numeric(metaclusters)) {
#     metaclusters <- list(metaclusters)
#   }
#
#   # Get all label coordinates
#   if (length(all_clust) == 1) {
#     mfis <- matrix(FlowSOM::GetClusterMFIs(input)[all_clust, channelpair],
#                   nrow = 1,
#                   dimnames = list(all_clust, channelpair))
#   } else {
#     mfis <- FlowSOM::GetClusterMFIs(input)[all_clust, channelpair]
#   }
#
#   df <- as.data.frame(mfis)
#
#   # Get label names
#   switch(label_type,
#          cluster = {labs <- rownames(df)},
#          metacluster = {labs <- input$metaclustering[all_clust]},
#          stop("label_type must be either 'cluster' or 'metacluster'")
#   )
#
#   print(df)
#
#   # Generate unlabeled ggplot2 plot list
#   p <- FlowSOM::Plot2DScatters(fsom = input,
#                               channelpairs = list(channelpair),
#                               clusters = clusters,
#                               metaclusters = metaclusters,
#                               plotFile = NULL,
#                               density = FALSE,
#                               centers = FALSE)
#
#   if (!is.null(metaclusters)) { # generates single plot  # color is normally 'black' but want to color label based on metacluster
#     p <- p[[length(p)]] + ggplot2::geom_label(data = df,
#                                               ggplot2::aes(x = df[,1], y = df[,2], label = labs),
#                                               color = labs, size = 2.5, fontface = "bold")
#   } else { # generates list of plots
#     for (i in 1:length(p)) {
#       p[[i]] <- p[[i]] + ggplot2::geom_label(data = df,
#                                              ggplot2::aes(x = df[,1], y = df[,2], label = labs),
#                                              color = "black", size = 2.5, fontface = "bold")
#     }
#   }
#
#   return(p)
# }

#' @keywords internal
#' @export
plotLabeled2DScatter.data.frame <- function(input, channelpair, clusters = NULL,
                                            metaclusters = NULL, labels = TRUE,
                                            colors = NULL) {
  Metacluster <- Cluster <- .data <- NULL
  only_clust <- TRUE

  if (is.null(clusters) && is.null(metaclusters)) {
    stop("clusters and/or metaclusters must be a vector of indices.")

  } else if (is.null(clusters) && !is.null(metaclusters)) {
    # Find all clusters within given metaclusters
    all_clust <- input %>%
      tidytable::filter(Metacluster %in% metaclusters) %>%
      tidytable::pull(Cluster) %>%
      unique()

    only_clust <- FALSE

  } else if (!is.null(clusters) && is.null(metaclusters)) {
    # Only clusters given
    all_clust <- clusters

  } else if (!is.null(clusters) && !is.null(metaclusters)) {
    # Combine nodes within given metaclusters with those specified by `clusters`
    meta_nodes <- input %>%
      tidytable::filter(Metacluster %in% metaclusters) %>%
      tidytable::pull(Cluster) %>%
      unique()
    all_clust <- unique(c(meta_nodes, clusters))

    only_clust <- FALSE
  }

  # Get all label coordinates
  if (length(all_clust) == 1) {
    clust_mfis <- as.matrix(getClusterMFIs(input, cols_to_use = channelpair)[all_clust, ],
                            nrow = 1,
                            dimnames = list(all_clust, channelpair))
  } else {
    clust_mfis <- getClusterMFIs(input, cols_to_use = channelpair)[all_clust, ]
  }

  # Create data frame for labels
  df_lab <- as.data.frame(clust_mfis)

  # Get the metacluster each cluster belongs to
  inds <- sapply(all_clust, function(clust) {match(clust, input$Cluster)})
  cluster_meta <- input %>%
    tidytable::slice(inds) %>%
    tidytable::mutate(all_clust) %>%
    tidytable::pull(Metacluster, name = all_clust)

  cluster_meta <- factor(cluster_meta)

  # Background data frame
  samp_ids <- sample.int(nrow(input), size = min(nrow(input), 5000))
  df_bg <- input %>%
    tidytable::select(channelpair) %>%
    tidytable::slice(samp_ids) %>%
    data.frame(check.names = FALSE)

  # Subset input
  if (only_clust) {
    # If only clusters given
    df_c <- input %>%
      tidytable::filter(Cluster %in% all_clust)
  } else {
    # If metaclusters given
    df_c <- input %>%
      tidytable::filter(Metacluster %in% metaclusters)
  }

  # Generate metacluster/cluster data frame of points
  df_c <- df_c %>%
    tidytable::select(tidytable::all_of(c(channelpair, "Metacluster"))) %>%
    tidytable::slice_sample(n = 3000, .by = "Metacluster") %>%
    data.frame(check.names = FALSE)

  labs <- rownames(df_lab)
  df_lab <- cbind(df_lab, Metacluster = cluster_meta)


  # Rename columns
  colnames(df_bg) <- colnames(df_c)[1:2] <- c("x1", "x2")

  # Generate plot
  p <- ggplot2::ggplot(data = df_c,
                       mapping = ggplot2::aes(x = .data$x1, y = .data$x2)) +
    ggplot2::geom_point(data = df_bg,
                        color = "grey",
                        size = 0.3) +
    ggplot2::geom_point(mapping = ggplot2::aes(color = .data$Metacluster),
                        size = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::xlab(channelpair[1]) +
    ggplot2::ylab(channelpair[2])

  # Add cluster labels if desired
  if (labels) {
    p <- p + ggplot2::geom_label(data = df_lab,
                                 mapping = ggplot2::aes(x = df_lab[,1], y = df_lab[,2], label = labs, color = .data$Metacluster),
                                 size = 2.5,
                                 fontface = "bold")
  }

  # Apply colors, if given
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }

  return(p)
}

#' plotUMAP
#'
#' Plots UMAP colored by metacluster.
#'
#' @param input A FlowSOM object or data table.
#' @param num_cells The number of cells to use for the dimension reduction.
#' Default is 5000.
#' @param labels Optional, a vector specifying the order in which metacluster
#' name will appear in the legend. All metacluster names present in \code{input}
#' should exists in this vector.
#' @param colors Optional, a vector of colors to use in the plot for each metacluster.
#' @param seed Optional, a seed for reproducibility.
#'
#' @return
#' A UMAP plot drawn with ggplot2.
#'
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
#' fsom <- FlowSOM::FlowSOM(file,
#'                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#'                          nClus = 10,
#'                          seed = 42,
#'                          xdim = 6,
#'                          ydim = 6)
#'
#' plotUMAP(fsom)
#'
#' plotUMAP(fsom)
plotUMAP <- function(input, num_cells = 5000, labels = NULL, colors = NULL, seed = NULL) {
  result <- UseMethod("plotUMAP")
  return(result)
}

#' @keywords internal
#' @export
plotUMAP.FlowSOM <- function(input, num_cells = 5000, labels = NULL,
                             colors = NULL, seed = NULL) {
  X1 <- X2 <- Metacluster <- NULL

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cols_used <- input$map$colsUsed
  meta_vec <- FlowSOM::GetMetaclusters(input)

  data <- input$data[, cols_used]

  nrows <- nrow(data)

  if (num_cells < nrows) {
    inds <- sample(nrows, num_cells)
  } else {
    inds <- 1:nrows
  }

  data <- data[inds, ]
  umap <- umap::umap(data)
  meta_vec <- meta_vec[inds]
  if (!is.null(labels)) {
    meta_vec <- factor(meta_vec, levels = labels)
  }

  umap_df <- data.frame(umap$layout, Metacluster = meta_vec, Indices = inds)

  # Draw plot
  p <- ggplot2::ggplot(umap_df) +
    scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2,
                                               color = Metacluster),
                                  pointsize = 2) +
    ggplot2::theme_void()

  if (!is.null(colors)) {
    colors <- colors
    names(colors) <- levels(input[["metaclustering"]])

    p <- p + ggplot2::scale_fill_manual(aesthetics = "color", values = colors)
  }

  return(p)
}

#' @keywords internal
#' @export
plotUMAP.data.frame <- function(input, num_cells = 5000, labels = NULL,
                                colors = NULL, seed = NULL) {
  X1 <- X2 <- Metacluster <- NULL

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cols_used <- attributes(input)$clustered
  meta_vec <- input$Metacluster

  input <- input %>%
    tidytable::select(cols_used)

  nrows <- nrow(input)

  if (num_cells < nrows) {
    inds <- sample(nrows, num_cells)
  } else {
    inds <- 1:nrows
  }

  data <- input[inds, ]
  umap <- umap::umap(data)
  meta_vec <- meta_vec[inds]
  if (!is.null(labels)) {
    meta_vec <- factor(meta_vec, levels = labels)
  }

  # add parameter for choosing order of legend

  umap_df <- data.frame(umap$layout, Metacluster = meta_vec, Indices = inds)

  # Draw plot
  p <- ggplot2::ggplot(umap_df) +
    scattermore::geom_scattermore(ggplot2::aes(x = X1, y = X2,
                                               color = Metacluster),
                                  pointsize = 2) +
    ggplot2::theme_void()

  if (!is.null(colors)) {
    colors <- colors
    names(colors) <- levels(input$Metacluster)

    p <- p + ggplot2::scale_fill_manual(aesthetics = "color", values = colors)
  }

  return(p)
}

#' plotTSNE
#'
#' Plot a t-SNE dimension reduction of FlowSOM data, colored by metacluster.
#'
#' @param fsom The FlowSOM object whose data you would like to plot.
#' @param num_cells The total number of cells you would like to use for the
#' dimension reduction.
#' @param point_size Optional parameter to adjust the size of the plotted points.
#' Default is 1.5.
#' @param seed Optional parameter to set a seed for reproducibility.
#'
#' @return A t-SNE plot made with ggplot2.
#' @export
#'
#' @examples
#' file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
#' fsom <- FlowSOM::FlowSOM(file,
#'                          colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
#'                          nClus = 10,
#'                          seed = 42,
#'                          xdim = 6,
#'                          ydim = 6)
#'
#' plotTSNE(fsom, num_cells = 2500)
#'
#' plotTSNE(fsom, num_cells = 2500, point_size = 3, seed = 42)
plotTSNE = function(fsom, num_cells, point_size = 1.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  tsne <- FlowSOM::PlotDimRed(fsom = fsom,
                             colorBy = "metaclusters",
                             colors = viridisLite::viridis(length(levels(fsom$metaclustering)), option = "turbo"),
                             cTotal = num_cells,
                             dimred = Rtsne::Rtsne)

  tsne$layers[[1]]$geom_params$pointsize <- point_size
  tsne$layers[[2]]$geom_params$max.overlaps <- 50

  return(tsne)
}

#' plotSampleProportions
#'
#' Plot cell type proportions by sample.
#'
#' @param input A data frame or table.
#'
#' @details
#' By default, the names used for each sample will be its file name.
#'
#' @return A stacked bar plot made with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
plotSampleProportions <- function(input) {
  rn <- Sample <- Proportion <- `Cell Type` <- NULL
    sample_prop <- makeCountMatrix(input,
                                 min_cells = 0,
                                 min_samples = 0) %>%
    tidytable::as_tidytable(.keep_rownames = TRUE) %>%
    tidytable::mutate(tidytable::across(tidytable::where(is.numeric), .fns = ~ ./sum(.)*100)) %>%
    tidytable::pivot_longer(cols = -rn,
                            names_to = "Sample",
                            values_to = "Proportion") %>%
    tidytable::rename(`Cell Type` = rn)

  ggplot2::ggplot(sample_prop,
                  ggplot2::aes(x = Sample, y = Proportion, fill = `Cell Type`)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Sample Names", y = "Proportion") +
    ggplot2::ggtitle("Cell Types by Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6)) +
    ggplot2::scale_fill_viridis_d(option = "turbo")
}

# reg_expr Optional, a regular expression to rename sample files.
# plotSampleProportionsOld = function(fsom, reg_expr = NULL) {
#   dir_prepr <- dir_prepr()
#   Proportion <- `Cell Type` <- NULL
#
#   sample_names <- list.files(dir_prepr)
#   sample_prcntgs <- c()
#   for (i in 1:length(sample_names)) {
#     print(paste0("Processing ", sample_names[i], "..."))
#     ind <- which(fsom$data[, "File"] == i)
#     fsom_subset <- FlowSOM::FlowSOMSubset(fsom, ind)
#     counts <- FlowSOM::GetCounts(fsom_subset, level = "metaclusters")
#     sample_prcntgs <- rbind(sample_prcntgs, counts/sum(counts))
#   }
#
#   if (!is.null(reg_expr)) {
#     sample_names <- sub(reg_expr, "\\1", sample_names)
#   }
#
#   sample_prcntgs <- cbind(sample_names, sample_prcntgs)
#   sample_prcntgs <- as.data.frame(sample_prcntgs)
#
#   prcntgs_long <- tidyr::pivot_longer(sample_prcntgs,
#                                       cols = -sample_names,
#                                       names_to = "Cell Type",
#                                       values_to = "Proportion")
#   prcntgs_long$Proportion <- methods::as(prcntgs_long$Proportion, "numeric")
#
#   ggplot2::ggplot(prcntgs_long,
#                   ggplot2::aes(x = sample_names, y = Proportion, fill = `Cell Type`)) +
#     ggplot2::geom_bar(stat = "identity") +
#     ggplot2::labs(x = "Sample Names", y = "Proportion") +
#     ggplot2::ggtitle("Cell Types by Sample") +
#     ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6)) +
#     ggplot2::scale_fill_viridis_d(option = "turbo")
# }

#' removeFilepath
#'
#' Helper to remove filepath from a filename.
#'
#' @param filenames A list of filenames.
#'
#' @return A list of filenames.
#'
#' @export
removeFilepath = function(filenames) {
  contains_path = grepl("/", filenames)
  new_filenames = ifelse(contains_path, sub(".*/", "", filenames), filenames)
  return(new_filenames)
}

#' plotClusterGroupProportions
#'
#' Plot group proportions by cluster in an MST or grid.
#'
#' @param fsom A FlowSOM object as generated by \code{\link[FlowSOM:FlowSOM]{FlowSOM()}}.
#' @param groups A list of lists, specifying the groups of interest and their file
#' names.
#' @param view A string specifying which version of the plot you would like returned
#' either "MST" or "grid".
#'
#' @details
#' For each cluster, a pie chart displaying the proportions
#' of each given group within said cluster is created. They
#' may be visualized as either an MST or a grid. The pie chart in
#' the bottom left may serve as a point of reference, as it shows what we would
#' expect the proportions to be if there was no difference between groups.
#'
#' @return An MST or grid plot.
#'
#' @export
plotClusterGroupProportions <- function(fsom, groups, view = c("MST", "grid")) {
  dir_clustr <- dir_clustr()

  file_names <- list.files(path = dir_clustr,
                           full.names = FALSE)
  all_group_files <- c()
  for (group in groups) {
    all_group_files <- c(all_group_files, group)
  }

  sample_names <- file_names[fsom$data[,"File"]]
  ind <- which(paste0(dir_clustr, sample_names) %in% all_group_files)

  if (length(ind) < length(sample_names)) {
    sample_names <- sample_names[ind]
    fsom <- FlowSOM::FlowSOMSubset(fsom, ind)
  }

  group_indicator <- rep("Unknown", length(sample_names))

  # check names for both "file_names" and "groups"
  for (i in 1:length(groups)) {
    groups[[i]] <- removeFilepath(groups[[i]])
    group_indicator[sample_names %in% unlist(groups[names(groups)[i]])] <- names(groups)[i]
  }

  file_pie_plot <- FlowSOM::PlotPies(fsom = fsom,
                                     cellTypes = factor(group_indicator),
                                     colorPalette = viridisLite::viridis(length(groups), option = "inferno"),
                                     equalNodeSize = TRUE,
                                     maxNodeSize = ifelse(view == "MST", 0.7, 1),
                                     view = view)

  group_proportions <- sapply(groups, length)
  group_proportions <- group_proportions / sum(group_proportions)

  angles <- cumsum(2 * pi * group_proportions)
  start_angles <- c(0, utils::head(angles, -1))

  file_pie_plot <- FlowSOM::AddStarsPies(p = file_pie_plot,
                                         arcs = data.frame(x0 = 0,
                                                           y0 = 0,
                                                           start = start_angles,
                                                           end = angles,
                                                           value = 1,
                                                           Markers = names(groups)),
                                         colorPalette = viridisLite::viridis(length(groups), option = "inferno")
  )

  return(file_pie_plot)
}

#' plotClusterFileProportions
#'
#' Plots file distributions by cluster in an MST or grid.
#'
#' @param fsom A FlowSOM object, as generated by \code{\link[FlowSOM:FlowSOM]{FlowSOM()}}.
#' @param view A string specifying which version of the plot you would like returned,
#' either "MST" or "grid".
#' @param reg_expr Optional, a regular expression to rename filenames.
#'
#' @details
#' The pie chart in the bottom left may serve as a point of reference,
#' as it shows what one would expect the proportions to be if
#' there was no difference between samples, and if the clustering
#' went well. Note that if a metacluster or cluster is made up
#' mostly of one sample's cells, this does not necessarily mean that the
#' clustering is flawed, and may instead indicate that the relevant sample has
#' unique biological characteristics.
#'
#' @return An MST or grid plot.
#'
#' @export
plotClusterFileProportions <- function(fsom, view = c("MST", "grid"), reg_expr = NULL) {
  sample_names <- list.files(path = dir_clustr(),
                             full.names = FALSE)

  if (!is.null(reg_expr)) {
    sample_names <- sub(reg_expr, "\\1", sample_names)
  }

  file_num <- length(sample_names)

  file_pie_plot <- FlowSOM::PlotPies(fsom = fsom,
                                     cellTypes = factor(sample_names[fsom$data[,"File"]]),
                                     colorPalette = viridisLite::viridis(file_num, option = "turbo"),
                                     equalNodeSize = TRUE,
                                     maxNodeSize = ifelse(view == "MST", 0.7, 1),
                                     view = view)

  proportions <- sapply(1:length(sample_names), function(i) {
    num_cells <- length(which(fsom$data[, "File"] == i))
    return(num_cells)
  })

  proportions <- proportions / sum(proportions)

  angles <- cumsum(2 * pi * proportions)
  start_angles <- c(0, utils::head(angles, -1))

  file_pie_plot <- FlowSOM::AddStarsPies(p = file_pie_plot,
                                         arcs = data.frame(x0 = 0,
                                                           y0 = 0,
                                                           start = start_angles,
                                                           end = angles,
                                                           value = 1,
                                                           Markers = sample_names),
                                         colorPalette = viridisLite::viridis(file_num, option = "turbo"))

  return(file_pie_plot)
}
