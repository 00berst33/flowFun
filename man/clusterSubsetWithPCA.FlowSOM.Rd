% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_functions.R
\name{clusterSubsetWithPCA.FlowSOM}
\alias{clusterSubsetWithPCA.FlowSOM}
\title{clusterSubsetWithPCA.FlowSOM}
\usage{
\method{clusterSubsetWithPCA}{FlowSOM}(
  input,
  pca_obj,
  num_components,
  num_clus,
  xdim = 10,
  ydim = 10,
  seed = NULL,
  convert_to_channels = TRUE,
  fsom_file = NULL
)
}
\description{
Recluster a subset of flow data using principal components of a filtered
aggregate .fcs file.
}
\keyword{internal}
