#' Print Method for BayesMallows Objects
#'
#' The default print method for a `BayesMallows` object.
#'
#' @param x An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param ... Other arguments passed to `print` (not used).
#'
#' @export
#'
#' @family posterior quantities
print.BayesMallows <- function(x, ...) {
  cat("Bayesian Mallows Model with", x$data$n_items, "items and", x$data$n_assessors, "assessors.\n")
  cat("Use functions assess_convergence() or plot() to visualize the object.")
}

#' @rdname print.BayesMallows
#' @export
print.BayesMallowsMixtures <- function(x, ...) {
  cat("Bayesian Mallows Mixtures Model with", length(x), "clusters.\n")
  cat("Use functions assess_convergence() or plot_elbow() to analyze.\n")
}

#' @rdname print.BayesMallows
#' @export
print.SMCMallows <- function(x, ...) {
  cat("Bayesian Mallows Model with", x$data$n_items, "items fitted with sequential Monte Carlo.\n")
  cat("Use the plot() to visualize the object.")
}
