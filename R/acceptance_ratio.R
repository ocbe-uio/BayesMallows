#' @title Get Acceptance Ratios
#' @description Extract acceptance ratio from Metropolis-Hastings
#'   algorithm used by [compute_mallows()] or by the move step in
#'   [update_mallows()] and [compute_mallows_sequentially()]. If
#'   burnin is not set in the call to [compute_mallows()], the acceptance ratio
#'   for all iterations will be reported. Otherwise the post burnin acceptance
#'   ratio is reported. For the SMC method the acceptance ratios apply to all
#'   iterations, since no burnin is needed in here.
#'
#' @param model_fit A model fit.
#' @param ... Other arguments passed on to other methods. Currently not used.
#'
#' @export
#' @return A list with elements `alpha_acceptance`, `rho_acceptance`, and
#'   `aug_acceptance`. Each element contains acceptance ratios (between 0 and 1)
#'   for the corresponding parameter proposals in the Metropolis-Hastings algorithm.
#'   For models with multiple chains, each element is a list with one acceptance
#'   ratio per chain. Higher values indicate higher acceptance rates for the
#'   Metropolis-Hastings proposals.
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
