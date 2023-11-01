#' Set prior parameters for Bayesian Mallows model
#'
#' @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the truncated exponential prior distribution of \eqn{\alpha}. Defaults
#'   to `0.1`. When `n_cluster > 1`, each mixture component
#'   \eqn{\alpha_{c}} has the same prior distribution.
#'
#' @param alpha_max Maximum value of `alpha` in the truncated exponential
#'   prior distribution.
#'
#' @param psi Positive integer specifying the concentration parameter \eqn{\psi}
#'   of the Dirichlet prior distribution used for the cluster probabilities
#'   \eqn{\tau_{1}, \tau_{2}, \dots, \tau_{C}}, where \eqn{C} is the value of
#'   `n_clusters`. Defaults to `10L`. When `n_clusters = 1`, this
#'   argument is not used.
#'
#' @param kappa_1 First shape parameter of the truncate Beta prior used for
#'   \eqn{theta} in the Bernoulli error model. Defaults to 1.0. See
#'   \insertCite{crispino2019}{BayesMallows} for details.
#'
#' @param kappa_2 Second shape parameter of the truncate Beta prior used for
#'   \eqn{theta} in the Bernoulli error model. Defaults to 1.0. See
#'   \insertCite{crispino2019}{BayesMallows} for details.
#'
#' @return An object of class `"BayesMallowsPriors"`, to be provided in the
#'   `priors` argument to [compute_mallows()].
#' @export
#'
#' @family modeling
#'
set_priors <- function(lambda = 0.001, alpha_max = 1e6, psi = 10,
                       kappa_1 = 1, kappa_2 = 1) {
  validate_positive(lambda)
  validate_positive(alpha_max)
  validate_integer(psi)
  validate_positive(psi)
  validate_positive(kappa_1)
  validate_positive(kappa_2)

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsPriors"
  ret
}
