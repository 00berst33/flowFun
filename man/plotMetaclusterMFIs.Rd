% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plotMetaclusterMFIs}
\alias{plotMetaclusterMFIs}
\title{plotMetaclusterMFIs}
\usage{
plotMetaclusterMFIs(input, cols_to_use = NULL, ...)
}
\arguments{
\item{input}{A FlowSOM object or data table.}

\item{cols_to_use}{A character vector specifying which markers you
would like to include in the plot. Default is the markers used for clustering.}

\item{...}{Additional parameters to pass to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}.}
}
\value{
A heatmap of MFIs made with \code{ComplexHeatmap::Heatmap()},
where each column is a marker of interest, and each row is a metacluster.
}
\description{
Plot a heatmap of metacluster MFIs.
}
\examples{
file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
fsom <- FlowSOM::FlowSOM(file,
                         colsToUse = c(10, 12:14, 16, 18:23, 25:32, 34),
                         nClus = 10,
                         seed = 42,
                         xdim = 6,
                         ydim = 6)

plotMetaclusterMFIs(fsom)
}
