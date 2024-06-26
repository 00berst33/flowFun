% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_functions.R
\name{getSampleMetaMatrix}
\alias{getSampleMetaMatrix}
\title{getSampleMetaMatrix}
\usage{
getSampleMetaMatrix(df_full, col_to_use)
}
\arguments{
\item{df_full}{A data frame with sample/metacluster information for the current
marker.}

\item{col_to_use}{Which column of the data frame should be used to create the
matrix.}
}
\value{
A matrix, where each column is a metacluster and each row is a sample.
The entries may be median expression on either a linear or logicle-transformed
scale.
}
\description{
Helper function for creating matrices of interest from \code{df_full}.
}
\keyword{internal}
