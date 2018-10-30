#' Print Method for BayesMallowsMixtures Objects
#'
#' The default print method for a \code{BayesMallowsMixtures} object.
#'
#'
#' @param x An object of type \code{BayesMallowsMixtures}, returned from
#'   \code{\link{compute_mallows_mixtures}}.
#'
#' @param ... Other arguments passed to \code{print} (not used).
#'
#' @export
#'
#'
print.BayesMallowsMixtures <- function(x, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because print.BayesMallowsMixtures must have the same
  # required arguments as base::print.

  if(!all(purrr::map_lgl(x, ~ inherits(.x, "BayesMallows")))) {
    stop("All elements of a BayesMallowsMixtures object must be of class BayesMallows.")
  }

  n_clusters <- purrr::map_int(x, ~ .x$n_clusters)

  cat("Collection of", length(x), "Bayesian Mallows Mixture Models with the following number of mixture components:\n",
      paste0(paste(n_clusters, collapse = ", "), "."))
}


