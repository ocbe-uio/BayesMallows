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
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
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
compute_mallows_mixtures <- function(n_clusters, ..., cl = NULL){
  stopifnot(is.numeric(n_clusters))
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  if(is.null(cl)){
    models <- purrr::map(n_clusters, function(x) {
      compute_mallows(..., n_clusters = x)
    })
  } else {
    args <- list(...)
    parallel::clusterExport(cl = cl, varlist = "args", envir = environment())
    models <- parallel::parLapply(cl = cl, X = n_clusters, fun = function(x){
      args$n_clusters <- x
      do.call(compute_mallows, args)
    })

  }

  class(models) <- "BayesMallowsMixtures"
  return(models)
}
