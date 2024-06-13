env1 <- new.env(parent = emptyenv())

# !!! edit setters !!!

env1$dir_rds <- file.path("RDS")
env1$dir_rds_unedited <- file.path("RDS", "Unedited")
env1$dir_rds_edited <- file.path("RDS", "Edited")
env1$dir_agg <- file.path("Data", "Aggregates")
env1$dir_result <- file.path("Analysis Results")
env1$dir_prepr <- file.path("Data", "Preprocessed Files")
env1$dir_clustr <- file.path("Data", "Clustered Files")
env1$dir_edger <- file.path("Analysis Results", "edgeR") 
env1$dir_limma <- file.path("Analysis Results", "limma")

#' Report dir_rds
#' @export
dir_rds <- function() {
  return(env1$dir_rds)
}

#' Edit dir_rds
#'
#' @param new_dir_rds A string defining the new directory for all RDS files.
#'
#' @export
set_dir_rds <- function(new_dir_rds) {
  old <- env1$dir_rds
  env1$dir_rds <- new_dir_rds
  invisible(old)
}

#' Report dir_rds_unedited
#' @export
dir_rds_unedited <- function() {
  return(env1$dir_rds_unedited)
}

#' Edit dir_rds_unedited
#'
#' @param new_dir_rds_unedited A string defining the new directory for unedited
#' RDS files.
#'
#' @export
set_dir_rds_unedited <- function(new_dir_rds_unedited) {
  old <- env1$dir_rds_unedited
  env1$dir_rds_unedited <- new_dir_rds_unedited
  invisible(old)
}

#' Report dir_rds_edited
#' @export
dir_rds_edited <- function() {
  return(env1$dir_rds_edited)
}

#' Edit dir_rds_edited
#'
#' @param new_dir_rds_edited A string defining the new directory for edited
#' RDS files.
#'
#' @export
set_dir_rds_edited <- function(new_dir_rds_edited) {
  old <- env1$dir_rds_edited
  env1$dir_rds_edited <- new_dir_rds_edited
  invisible(old)
}

#' Report dir_agg
#' @export
dir_agg <- function() {
  return(env1$dir_agg)
}

#' Edit dir_agg
#'
#' @param new_dir_agg A string defining the new directory for aggregate files.
#'
#' @export
set_dir_agg <- function(new_dir_agg) {
  old <- env1$dir_agg
  env1$dir_agg <- new_dir_agg
  invisible(old)
}

#' Report dir_result
#' @export
dir_result <- function() {
  return(env1$dir_result)
}

#' Edit dir_result
#'
#' @param new_dir_result A string defining the new directory for analysis results.
#'
#' @export
set_dir_result <- function(new_dir_result) {
  old <- env1$dir_result
  env1$dir_result <- new_dir_result
  invisible(old)
}

#' Report dir_prepr
#' @export
dir_prepr <- function() {
  return(env1$dir_prepr)
}

#' Edit dir_prepr
#'
#' @param new_dir_prepr A string defining the new directory for preprocessed .fcs
#' files.
#'
#' @export
set_dir_prepr <- function(new_dir_prepr) {
  old <- env1$dir_prepr
  env1$dir_prepr <- new_dir_prepr
  invisible(old)
}

#' Report dir_clustr
#' @export
dir_clustr <- function() {
  return(env1$dir_clustr)
}

#' Edit dir_clustr
#'
#' @param new_dir_clustr A string defining the directory to store results of
#' differential abundance analysis from edgeR.
#'
#' @export
set_dir_clustr <- function(new_dir_clustr) {
  old <- env1$dir_clustr
  env1$dir_clustr <- new_dir_clustr
  invisible(old)
}

#' Report dir_edger
#' @export
dir_edger <- function() {
  return(env1$dir_edger)
}

#' Edit dir_edger
#'
#' @param new_dir_edger A string defining the directory to store results of
#' differential abundance analysis from edgeR.
#'
#' @export
set_dir_edger <- function(new_dir_edger) {
  old <- env1$dir_edger
  env1$dir_edger <- new_dir_edger
  invisible(old)
}

#' Report dir_limma
#' @export
dir_limma <- function() {
  return(env1$dir_limma)
}

#' Edit dir_limma
#'
#' @param new_dir_limma A string defining the directory to store results of
#' differential expression analysis from limma.
#'
#' @export
set_dir_limma <- function(new_dir_limma) {
  old <- env1$dir_limma
  env1$dir_limma <- new_dir_limma
  invisible(old)
}