% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.R
\name{editFlowSOMMetaData}
\alias{editFlowSOMMetaData}
\title{editFlowSOMMetaData}
\usage{
editFlowSOMMetaData(fsom, metadata, clear = NULL)
}
\arguments{
\item{fsom}{A FlowSOM object to annotate.}

\item{metadata}{A named list of metadata to add to the FlowSOM object.}

\item{clear}{Boolean, should the existing metadata be cleared before editing.}
}
\value{
A FlowSOM object with edited metadata.
}
\description{
editFlowSOMMetaData
}
\details{
It is useful to keep track of:
\itemize{
\item aggregate file used for clustering
\itemize{
\item number of cells aggregated
\item which files were used and in what order they were read
}
\item which directory contains the clustered files resulting from the FlowSOM object
\item any parent FlowSOM objects, if current FlowSOM object is a reclustering
}
}
\keyword{internal}
