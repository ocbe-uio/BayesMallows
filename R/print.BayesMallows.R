#' Print Method for BayesMallows Objects
#'
#' The default print method for a \code{BayesMallows} object.
#'
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param ... Other arguments passed to \code{print} (not used).
#'
#' @export
#'
#'
print.BayesMallows <- function(x, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because print.BayesMallows must have the same
  # required arguments as base::print.

  if(is.null(x$n_items) || is.null(x$n_assessors)) {
    stop("BayesMallows object must have elements n_items and n_assessors.")
  }
  cat("Bayesian Mallows Model with", x$n_items, "items and", x$n_assessors, "assessors.\n")
  cat("Use functions assess_convergence() or plot() to visualize the object.")
}


