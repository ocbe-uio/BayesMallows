#' Sample from prior distribution
#'
#' Function to obtain samples from the prior distributions of the Bayesian
#' Mallows model. Intended to be given to [update_mallows()].
#'
#' @param n An integer specifying the number of samples to take.
#' @param n_items An integer specifying the number of items to be ranked.
#' @param priors An object of class "BayesMallowsPriors" returned from
#'   [set_priors()].
#'
#' @return An object of class "BayesMallowsPriorSample", containing `n`
#' independent samples of \eqn{\alpha} and \eqn{\rho}.
#'
#' @export
#'
#' @family modeling
#' @example /inst/examples/sample_prior_example.R
sample_prior <- function(n, n_items, priors = set_priors()) {
  validate_positive(n)
  validate_positive(n_items)
  ret <- list(
    alpha = rexp(n, rate = priors$lambda),
    rho = replicate(n, sample(n_items, n_items)),
    priors = priors,
    n_items = n_items,
    items = seq_len(n_items)
  )
  class(ret) <- "BayesMallowsPriorSamples"
  ret
}
