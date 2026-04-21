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
    if (!is.null(attr(input, "clustered"))) {
      cols_to_use <- attr(input, "clustered")
    } else {
      stop("Default not found, please specify `cols_to_use`")
    }
  }
  # If input is only channel name
  if (is.character(cols_to_use) & !any(cols_to_use %in% colnames(input))) {
    marker_cols <- sub(".*<(.*)>.*", "\\1", colnames(input))
    match_idx <- match(cols_to_use, marker_cols)
    # Set channelpair to pretty column names
    cols_to_use <- colnames(input)[match_idx]
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
    if (!is.null(attr(input, "clustered"))) {
      cols_to_use <- attr(input, "clustered")
    } else {
      stop("Default not found, please specify `cols_to_use`")
    }
  }
  # If input is only channel name
  if (is.character(cols_to_use) & !any(cols_to_use %in% colnames(input))) {
    marker_cols <- sub(".*<(.*)>.*", "\\1", colnames(input))
    match_idx <- match(cols_to_use, marker_cols)
    # Set channelpair to pretty column names
    cols_to_use <- colnames(input)[match_idx]
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
annotateMFIHeatmap.data.frame <- function(merged_input, original_input = NULL, cols_to_use = NULL, ...) {
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
#' # not run
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

  # If channelpair not found in table column names
  if (!any(channelpair %in% colnames(input))) {
    marker_cols <- sub(".*<(.*)>.*", "\\1", colnames(input))
    match_idx <- match(channelpair, marker_cols)
    # Set channelpair to pretty column names
    channelpair <- colnames(input)[match_idx]
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


#' @keywords internal
#' @export
plotLabeled2DScatter.FlowSOM <- function(input, channelpair, clusters = NULL, # needs edit for one case
                                         metaclusters = NULL,
                                         label_type = c("cluster", "metacluster")) {

  if (is.null(clusters) && is.null(metaclusters)) {
    stop("clusters and/or metaclusters must be a vector of indices.")

  } else if (is.null(clusters) && !is.null(metaclusters)) {
    # color by metacluster, plot colored labels with cluster numbers
    all_clust <- which(input$metaclustering %in% metaclusters)

  } else if (!is.null(clusters) && is.null(metaclusters)) {
    # color by cluster
    all_clust <- clusters

  } else if (!is.null(clusters) && !is.null(metaclusters)) {
    # color by metacluster, plot colored labels with cluster numbers
    meta_nodes <- which(input$metaclustering %in% metaclusters)
    all_clust <- unique(c(meta_nodes, clusters))

  }

  # edit !!!
  # values in `metaclusters` and levels(input$metaclustering) can conflict
  if (is.character(metaclusters)) {
    metaclusters <- which(levels(input$metaclustering) %in% metaclusters)
  } else if (is.numeric(metaclusters)) {
    metaclusters <- list(metaclusters)
  }

  # Get all label coordinates
  if (length(all_clust) == 1) {
    mfis <- matrix(FlowSOM::GetClusterMFIs(input)[all_clust, channelpair],
                  nrow = 1,
                  dimnames = list(all_clust, channelpair))
  } else {
    mfis <- FlowSOM::GetClusterMFIs(input)[all_clust, channelpair]
  }

  df <- as.data.frame(mfis)

  # Get label names
  switch(label_type,
         cluster = {labs <- rownames(df)},
         metacluster = {labs <- input$metaclustering[all_clust]},
         stop("label_type must be either 'cluster' or 'metacluster'")
  )

  print(df)

  # Generate unlabeled ggplot2 plot list
  p <- FlowSOM::Plot2DScatters(fsom = input,
                              channelpairs = list(channelpair),
                              clusters = clusters,
                              metaclusters = metaclusters,
                              plotFile = NULL,
                              density = FALSE,
                              centers = FALSE)

  if (!is.null(metaclusters)) { # generates single plot  # color is normally 'black' but want to color label based on metacluster
    p <- p[[length(p)]] + ggplot2::geom_label(data = df,
                                              ggplot2::aes(x = df[,1], y = df[,2], label = labs),
                                              color = labs, size = 2.5, fontface = "bold")
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


#' plotSampleProportions
#'
#' Plot cell type proportions by sample.
#'
#' @param count_mat A data frame or table, where rows are clusters and
#' columns are samples.
#'
#' @details
#' By default, the names used for each sample will be its file name.
#'
#' @return A stacked bar plot made with \code{\link[ggplot2]{ggplot2}}.
#' @seealso [getExprMatDE()]
#'
#' @export
plotSampleProportions <- function(count_mat) {
  rn <- Sample <- Proportion <- `Cell Type` <- NULL
    sample_prop <- count_mat %>%
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

#' plotMarkerDensityByChannel
#'
#' @param gs A `GatingSet`
#' @param channel The channel to plot densities for
#' @param population A `character` indicating which subpopulation to use
#' @param inverse A boolean, whether or not the data should be inverse transformed
#' before plotting. `FALSE` by default
#' @param verbose Boolean specifying whether or not to print progress updates as
#' function runs. Default is `TRUE`
#'
#' @return A faceted plot with density plots of the chosen channel for each sample
#' @export
plot1DMarkerDensities <- function(gs, channel, population = "root",
                                  facet_by = c("samples", "subpopulations"),
                                  inverse = FALSE, verbose = TRUE) {
  channel <- rlang::enquo(channel)

  if (facet_by == "samples") {
    # Plot marker expression by sample
    p <- ggcyto::ggcyto(gs, ggplot2::aes(x = !!channel), subset = population) +
      ggplot2::geom_density()
    if (inverse) {
      p <- p + ggcyto::axis_x_inverse_trans()
    }
  }

  if (facet_by == "subpopulations") {
    # Plot by metacluster
    meta_pops <- flowWorkspace::gs_pop_get_children(gs, y = population, path = "auto")
    # Make plots for each metacluster
    p <- lapply(meta_pops, function(pop) {
      str <- paste0("Plotting marker density for metacluster named ", pop, "...")
      print(str)
      p <- ggcyto::ggcyto(gs, ggplot2::aes(x = !!channel), subset = pop) +
        ggplot2::geom_density() +
        ggplot2::facet_null()
      if (inverse) {
        p <- p + ggcyto::axis_x_inverse_trans()
      }
      p <- ggcyto::as.ggplot(p)
      return(p)
    })
    p <- patchwork::wrap_plots(unlist(p))
  }

  return (p)
}
