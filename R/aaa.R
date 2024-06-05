env1 <- new.env(parent = emptyenv())

env1$dir_rds = "RDS/"
env1$dir_rds_unedited = "RDS/Unedited/"
env1$dir_rds_edited = "RDS/Edited/"
env1$dir_agg = "Data/Aggregates/"
env1$dir_result = "Analysis Results/"
env1$dir_prepr = "Data/Preprocessed Raw/"
env1$dir_clustr = "Data/Clustered Raw/"

#' Report dir_rds
#' @export
dir_rds <- function() {
  return(env1$dir_rds)
}

#' Edit dir_rds
#' @export
set_dir_rds <- function(new_dir_rds) {
  assign("dir_rds", new_dir_rds, env = env1)
}

#' Report dir_rds_unedited
#' @export
dir_rds_unedited <- function() {
  return(env1$dir_rds_unedited)
}

#' Edit dir_rds_unedited
#' @export
set_dir_rds_unedited <- function(new_dir_rds_unedited) {
  assign("dir_rds_unedited", new_dir_rds_unedited, env = env1)
}

#' Report dir_rds_edited
#' @export
dir_rds_edited <- function() {
  return(env1$dir_rds_edited)
}

#' Edit dir_rds_edited
#' @export
set_dir_rds_edited <- function(new_dir_rds_edited) {
  assign("dir_rds_edited", new_dir_rds_edited, env = env1)
}

#' Report dir_agg
#' @export
dir_agg <- function() {
  return(env1$dir_agg)
}

#' Edit dir_agg
#' @export
set_dir_agg <- function(new_dir_agg) {
  assign("dir_agg", new_dir_agg, env = env1)
}

#' Report dir_result
#' @export
dir_result <- function() {
  return(env1$dir_result)
}

#' Edit dir_result
#' @export
set_dir_result <- function(new_dir_result) {
  assign("dir_result", new_dir_result, env = env1)
}

#' Report dir_prepr
#' @export
dir_prepr <- function() {
  return(env1$dir_prepr)
}

#' Edit dir_prepr
#' @export
set_dir_prepr <- function(new_dir_prepr) {
  assign("dir_prepr", new_dir_prepr, env = env1)
}

#' Report dir_clustr
#' @export
dir_clustr <- function() {
  return(env1$dir_clustr)
}

#' Edit dir_clustr
#' @export
set_dir_clustr <- function(new_dir_clustr) {
  assign("dir_clustr", new_dir_clustr, env = env1)
}
