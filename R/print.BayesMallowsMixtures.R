#' Print Method for BayesMallowsMixtures Objects
#'
#' The default print method for a `BayesMallowsMixtures` object.
#'
#'
#' @param x An object of type `BayesMallowsMixtures`, returned from
#'   [compute_mallows_mixtures()].
#'
#' @param ... Other arguments passed to `print` (not used).
#'
#' @export
#'
#' @family posterior quantities
print.BayesMallowsMixtures <- function(x, ...) {
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because print.BayesMallowsMixtures must have the same
  # required arguments as base::print.

  if (!Reduce(`&`, lapply(x, function(x) inherits(x, "BayesMallows")))) {
    stop("All elements of a BayesMallowsMixtures object must be of class BayesMallows.")
  }

  n_clusters <- vapply(x, function(x) x$n_clusters, integer(1))

  cat(
    "Collection of", length(x), "Bayesian Mallows Mixture Models with the following number of mixture components:\n",
    paste0(paste(n_clusters, collapse = ", "), ".")
  )
}
