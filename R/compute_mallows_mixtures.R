#' Compute Mixtures of Mallows Models
#'
#' Convenience function for computing Mallows models with varying
#' numbers of mixtures. This is useful for deciding the number of
#' mixtures to use in the final model.
#'
#' @param n_clusters Integer vector specifying the number of clusters
#' to use.
#'
#' @param ... Other named arguments, passed to \code{\link{compute_mallows}}.
#'
#' @return A list of Mallows models of class \code{BayesMallowsMixtures}, with one element
#' for each number of mixtures that
#' was computed. This object can be studied with \code{\link{plot_elbow}}.
#'
#'
#' @export
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#'
compute_mallows_mixtures <- function(n_clusters, ...){
  stopifnot(is.numeric(n_clusters))

  models <- purrr::map(n_clusters, function(x) {
    message(paste0("Computing Mallows model with ", x, " clusters."))
    compute_mallows(..., n_clusters = x)
  })

  class(models) <- "BayesMallowsMixtures"
  return(models)
}
