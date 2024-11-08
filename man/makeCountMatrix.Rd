% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_functions.R
\name{makeCountMatrix}
\alias{makeCountMatrix}
\title{makeCountMatrix}
\usage{
makeCountMatrix(
  input,
  meta_names = NULL,
  min_cells = 3,
  min_samples = NULL,
  ...
)
}
\arguments{
\item{input}{A FlowSOM object or data frame.}

\item{meta_names}{A vector of metacluster names of interest. By default, all
metacluster names are used.}

\item{min_cells}{An integer, the minimum number of cells a metacluster should
have in a specified number of samples to be included in the analysis.}

\item{min_samples}{An integer, the minimum number of samples a metacluster
should have at least \code{min_cells} events in to be included in the analysis.
By default, this is half the total number of samples.}

\item{...}{other arguments:

\code{sample_df}: If \code{input} is a FlowSOM object, a data frame from
\code{\link[=prepareSampleInfo]{prepareSampleInfo()}}.}
}
\value{
A matrix, where each column represents a sample, and each row
represents a metacluster.
}
\description{
Generate matrix of sample/metacluster cell counts.
}
\seealso{
\code{\link[=makeCountMatrix.FlowSOM]{makeCountMatrix.FlowSOM()}}
}
