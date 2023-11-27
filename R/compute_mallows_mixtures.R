#' Compute Mixtures of Mallows Models
#'
#' Convenience function for computing Mallows models with varying numbers of
#' mixtures. This is useful for deciding the number of mixtures to use in the
#' final model.
#'
#' @param n_clusters Integer vector specifying the number of clusters to use.
#' @inheritParams compute_mallows
#'
#'
#' @return A list of Mallows models of class `BayesMallowsMixtures`, with
#'   one element for each number of mixtures that was computed. This object can
#'   be studied with [plot_elbow()].
#'
#' @details
#' The `n_clusters` argument to [set_model_options()] is ignored
#' when calling `compute_mallows_mixtures`.
#'
#'
#' @family modeling
#'
#' @export
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#'
compute_mallows_mixtures <- function(
    n_clusters,
    data,
    model_options = set_model_options(),
    compute_options = set_compute_options(),
    priors = set_priors(),
    initial_values = set_initial_values(),
    logz_estimate = NULL,
    verbose = FALSE,
    cl = NULL) {
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  if (is.null(cl)) {
    lapplyfun <- lapply
  } else {
    varlist <- c(
      "data", "model_options", "compute_options", "priors", "initial_values",
      "logz_estimate", "verbose"
    )

    parallel::clusterExport(cl = cl, varlist = varlist, envir = environment())
    lapplyfun <- function(X, FUN, ...) {
      parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    }
  }

  models <- lapplyfun(n_clusters, function(x) {
    model_options$n_clusters <- x
    compute_mallows(
      data = data, model_options = model_options, compute_options = compute_options,
      priors = priors, initial_values = initial_values, verbose = verbose
    )
  })

  class(models) <- "BayesMallowsMixtures"
  return(models)
}
