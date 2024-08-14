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


#' plotClusterMFIs
#'
#' Plot a heatmap of cluster MFIs.
#'
#' @param fsom The FlowSOM object whose data you would like to plot.
#' @param markers_of_interest A character vector specifying which markers you
#' would like to include in the plot. Default is \code{fsom$map$colsUsed} (the markers
#' used for clustering).
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
plotClusterMFIs = function(fsom, markers_of_interest = fsom$map$colsUsed,
                           metaclusters = NULL, ...) {

  # Get cluster MFIs for each marker/channel of interest
  mfi_mat <- FlowSOM::GetClusterMFIs(fsom, colsUsed = FALSE, prettyColnames = FALSE)
  mfi_mat <- mfi_mat[, which(colnames(mfi_mat) %in% FlowSOM::GetChannels(fsom, markers_of_interest))]

  # Plot only clusters belonging to given metaclusters
  if (!is.null(metaclusters)) {
    ind <- fsom$metaclustering %in% metaclusters
    if (!(TRUE %in% ind)) {
      stop("No clusters found in the given metacluster(s). Check that either your indices are within bounds, or names have no typos.")
    }
    mfi_mat <- mfi_mat[ind, ]
  }

  colnames(mfi_mat) <- fsom$prettyColnames[colnames(mfi_mat)]

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
#' Plot 2D scatterplots of FlowSOM clusters, with (meta)cluster number labels for
#' each cluster center.
#'
#' @param fsom A FlowSOM object as generated by the FlowSOM function.
#' @param channelpair A vector of two channel names.
#' @param clusters A vector specifying indices of clusters of interest.
#' @param metaclusters A vector specifying indices of metaclusters of interest.
#' @param label_type A string specifying whether you would like the cluster
#' centers to be labeled with the cluster or metacluster they belong to, either
#' "cluster" or "metacluster".
#'
#' @details
#' It is only necessary to define one of the arguments \code{clusters} and
#' \code{metaclusters}, but the user may use both if they wish.
#'
#' Note that the output of this function changes slightly depending on the combination
#' of clusters and metaclusters the user is plotting. If the user plots metaclusters,
#' with or without additional clusters defined in the \code{clusters} argument,
#' then the output will be a single ggplot object, where each metacluster has a
#' unique color. If the user plots only clusters, such that the \code{metaclusters}
#' argument remains \code{NULL}, then a list of ggplot objects will be returned.
#' Each plot will correspond to a cluster, and color only the cells belonging to
#' said cluster. See examples for further explanation.
#'
#' The (meta)cluster number labels are included on every plot.
#'
#' @return A list of scatterplots or single scatterplot made with ggplot and
#' [FlowSOM::Plot2DScatters()].
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
#' # Returns a single plot colored by metacluster
#' plotLabeled2DScatter(fsom,
#'                      channelpair = c("APC-Cy7-A", "BUV615-P-A"),
#'                      metaclusters = c(1, 2, 7),
#'                      label_type = "metacluster")
#'
#' # Also returns a single plot colored by metacluster
#' plotLabeled2DScatter(fsom,
#'                      channelpair = c("APC-Cy7-A", "BUV805-A"),
#'                      clusters = c(25, 33),
#'                      metaclusters = c(1, 4),
#'                      label_type = "cluster")
#'
#' # Returns a list of two plots, once for each cluster
#' p = plotLabeled2DScatter(fsom,
#'                          channelpair = c("APC-Cy7-A", "BUV805-A"),
#'                          clusters = c(4, 19),
#'                          label_type = "cluster")
#' print(p[[1]])
#' print(p[[2]])
plotLabeled2DScatter = function(fsom, channelpair, clusters = NULL, metaclusters = NULL,
                                label_type = c("cluster", "metacluster")) {

  if (is.null(clusters) && is.null(metaclusters)) {
    stop("clusters and/or metaclusters must be a vector of indices.")

  } else if (is.null(clusters) && !is.null(metaclusters)) {
    all_clust <- which(fsom$metaclustering %in% metaclusters)

  } else if (!is.null(clusters) && is.null(metaclusters)) {
    all_clust <- clusters

  } else if (!is.null(clusters) && !is.null(metaclusters)) {
    meta_nodes <- which(fsom$metaclustering %in% metaclusters)
    all_clust <- unique(c(meta_nodes, clusters))

  }

  # edit !!!
  # values in `metaclusters` and levels(fsom$metaclustering) can conflict
  if (is.character(metaclusters)) {
    metaclusters <- which(levels(fsom$metaclustering) %in% metaclusters)
  } else if (is.numeric(metaclusters)) {
    metaclusters <- list(metaclusters)
  }

  # Get all label coordinates
  if (length(all_clust) == 1) {
    mfis <- matrix(FlowSOM::GetClusterMFIs(fsom)[all_clust, channelpair],
                  nrow = 1,
                  dimnames = list(all_clust, channelpair))
  } else {
    mfis <- FlowSOM::GetClusterMFIs(fsom)[all_clust, channelpair]
  }

  df <- as.data.frame(mfis)

  # Get label names
  switch(label_type,
         cluster = {labs <- rownames(df)},
         metacluster = {labs <- fsom$metaclustering[all_clust]},
         stop("label_type must be either 'cluster' or 'metacluster'")
  )

  print(metaclusters)

  # Generate unlabeled ggplot2 plot list
  p <- FlowSOM::Plot2DScatters(fsom = fsom,
                              channelpairs = list(channelpair),
                              clusters = clusters,
                              metaclusters = metaclusters,
                              plotFile = NULL,
                              density = FALSE,
                              centers = FALSE)

  if (!is.null(metaclusters)) { # generates single plot
    p <- p[[length(p)]] + ggplot2::geom_label(data = df,
                                             ggplot2::aes(x = df[,1], y = df[,2], label = labs),
                                             color = "black", size = 2.5, fontface = "bold")
  } else { # generates list of plots
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + ggplot2::geom_label(data = df,
                                            ggplot2::aes(x = df[,1], y = df[,2], label = labs),
                                            color = "black", size = 2.5, fontface = "bold")
    }
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

#' plotUMAP.FlowSOM
#'
#' *parameter tweaking
#'
#' @param input A data table or FlowSOM object.
#' @param num_cells The number of cells to use to generate the UMAP. Default is 5000.
#' @param seed Optional, a seed for reproducibility.
#'
#' @return A UMAP plot drawn with ggplot2.
#'
#' @export
plotUMAP.FlowSOM <- function(input, num_cells = 5000, labels = NULL,
                             colors = NULL, seed = NULL) {
  X1 <- X2 <- Metacluster <- NULL

  if (!is.null(seed)) {
    set.seed(seed)
  }

  cols_used <- input$map$colsUsed
  meta_vec <- FlowSOM::GetMetaclusters(input)

  input <- input$data[, cols_used]

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
    names(colors) <- levels(input$metaclustering)

    p <- p + ggplot2::scale_fill_manual(aesthetics = "color", values = colors)
  }

  return(p)
}

#' plotUMAP.data.frame
#'
#' @param input A data table or FlowSOM object.
#' @param num_cells The number of cells to use to generate the UMAP. Default is 5000.
#' @param seed Optional, a seed for reproducibility.
#'
#' @return A UMAP plot drawn with ggplot2.
#'
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
#' @param fsom A FlowSOM object as generated by \code{\link[FlowSOM:FlowSOM]{FlowSOM()}},
#' [clusterSubset()], or [clusterSubsetWithPCA()].
#' @param reg_expr Optional, a regular expression to rename sample files.
#'
#' @details
#' By default, the names used for each sample will be its file name. However,
#' you may rename the samples by providing a regular expression with a capturing
#' group, which will be applied to each file name.
#'
#' @return A stacked bar plot made with \code{\link[ggplot2]{ggplot2}}.
#'
#' @export
plotSampleProportions = function(fsom, reg_expr = NULL) {
  dir_prepr <- dir_prepr()
  Proportion <- `Cell Type` <- NULL

  sample_names <- list.files(dir_prepr)
  sample_prcntgs <- c()
  for (i in 1:length(sample_names)) {
    print(paste0("Processing ", sample_names[i], "..."))
    ind <- which(fsom$data[, "File"] == i)
    fsom_subset <- FlowSOM::FlowSOMSubset(fsom, ind)
    counts <- FlowSOM::GetCounts(fsom_subset, level = "metaclusters")
    sample_prcntgs <- rbind(sample_prcntgs, counts/sum(counts))
  }

  if (!is.null(reg_expr)) {
    sample_names <- sub(reg_expr, "\\1", sample_names)
  }

  sample_prcntgs <- cbind(sample_names, sample_prcntgs)
  sample_prcntgs <- as.data.frame(sample_prcntgs)

  prcntgs_long <- tidyr::pivot_longer(sample_prcntgs,
                                      cols = -sample_names,
                                      names_to = "Cell Type",
                                      values_to = "Proportion")
  prcntgs_long$Proportion <- methods::as(prcntgs_long$Proportion, "numeric")

  ggplot2::ggplot(prcntgs_long,
                  ggplot2::aes(x = sample_names, y = Proportion, fill = `Cell Type`)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::labs(x = "Sample Names", y = "Proportion") +
    ggplot2::ggtitle("Cell Types by Sample") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 6)) +
    ggplot2::scale_fill_viridis_d(option = "turbo")
}

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
