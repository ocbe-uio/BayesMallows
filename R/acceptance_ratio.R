#' @title Get Acceptance Ratios
#' @description Extract acceptance ratio from Metropolis-Hastings
#'   algorithm used by [compute_mallows()] or by the move step in
#'   [update_mallows()] and [compute_mallows_sequentially()]. Currently the
#'   function only returns the values, but it will be refined in the future. If
#'   burnin is not set in the call to [compute_mallows()], the acceptance ratio
#'   for all iterations will be reported. Otherwise the post burnin acceptance
#'   ratio is reported. For the SMC method the acceptance ratios apply to all
#'   iterations, since no burnin is needed in here.
#'
#' @param model_fit A model fit.
#' @param ... Other arguments passed on to other methods. Currently not used.
#'
#' @export
#' @example /inst/examples/get_acceptance_ratios_example.R
#'
#' @family posterior quantities
#'
get_acceptance_ratios <- function(model_fit, ...) {
  UseMethod("get_acceptance_ratios")
}

#' @export
#' @rdname get_acceptance_ratios
get_acceptance_ratios.BayesMallows <- function(model_fit, ...) {
  model_fit$acceptance_ratios
}

#' @export
#' @rdname get_acceptance_ratios
get_acceptance_ratios.SMCMallows <- function(model_fit, ...) {
  model_fit$acceptance_ratios
}
