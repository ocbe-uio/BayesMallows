.onUnload <- function(libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @noRd
validate_permutation <- function(vec) {
  if (!any(is.na(vec))) {
    return(all(sort(vec) == seq_along(vec)))
  } else if (all(is.na(vec))) {
    return(TRUE)
  } else {
    return(all(vec[!is.na(vec)] <= length(vec)) &&
      all(vec[!is.na(vec)] >= 1) && !any(duplicated(vec[!is.na(vec)])))
  }
}

count_jobs_per_cluster <- function(n_iterations, n_clusters) {
  # Split n_iterations into each cluster
  n_iterations_vec <- rep(floor(n_iterations / n_clusters), n_clusters)
  i <- 1
  while (sum(n_iterations_vec) != n_iterations) {
    n_iterations_vec[i] <- n_iterations_vec[i] + 1
    if (i > n_clusters) break
  }
  n_iterations_vec
}
