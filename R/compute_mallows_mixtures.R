#' Compute Mixtures of Mallows Models
#'
#' Convenience function for computing Mallows models with varying
#' numbers of mixtures. This is useful for determing the number of
#' mixtures to use in the final model.
#'
#' @param n_clusters Integer vector specifying the number of clusters
#' to use.
#'
#' @param ... Other named arguments, passed to \code{\link{compute_mallows}}.
#'
#' @return A list of Mallows model, one for each number of mixtures that
#' was computed. This object can be studied with \code{\link{plot_elbow}}.
#'
#' @details \code{compute_mallows_mixtures} always sets \code{include_wcd = TRUE}
#' in the code to \code{\link{compute_mallows}}, since the purpose of this
#' function is to create an elbow plot with \code{\link{plot_elbow}}.
#'
#' @export
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#'
compute_mallows_mixtures <- function(n_clusters, ...){
  stopifnot(is.numeric(n_clusters))

  models <- purrr::map(n_clusters, function(x) {
    compute_mallows(..., n_clusters = x, include_wcd = TRUE)
  })
}
