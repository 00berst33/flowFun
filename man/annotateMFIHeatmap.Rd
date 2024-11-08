% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{annotateMFIHeatmap}
\alias{annotateMFIHeatmap}
\title{annotateMFIHeatmap}
\usage{
annotateMFIHeatmap(
  merged_input,
  original_input = NULL,
  cols_to_use = NULL,
  ...
)
}
\arguments{
\item{merged_input}{The final table or FlowSOM object, whose metaclusters have
been merged and named via either \code{\link[=editTableMetaclusters]{editTableMetaclusters()}} or
\code{\link[FlowSOM:UpdateMetaclusters]{FlowSOM::UpdateMetaclusters()}}.}

\item{original_input}{If \code{merged_input} is a FlowSOM object, the initial
FlowSOM object before merging should be provided here.}

\item{cols_to_use}{A character vector specifying which markers you
would like to include in the plot. Default is the markers used for clustering.}

\item{...}{Additional parameters to pass to \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}.}
}
\value{
An annotated heatmap generated by \code{\link[ComplexHeatmap:Heatmap]{ComplexHeatmap::Heatmap()}}.
}
\description{
Annotate an unmerged heatmap of MFIs with each metacluster's final assignment.
}
\details{
Additional details...
}
