#' @title Set prior parameters for Bayesian Mallows model
#'
#' @description Set values related to the prior distributions for the Bayesian
#'   Mallows model.
#'
#' @param gamma Strictly positive numeric value specifying the shape parameter
#'   of the gamma prior distribution of \eqn{\alpha}. Defaults to `1`, thus
#'   recovering the exponential prior distribution used by
#'   \insertCite{vitelli2018}{BayesMallows}.
#'
#' @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the gamma prior distribution of \eqn{\alpha}. Defaults
#'   to `0.001`. When `n_cluster > 1`, each mixture component \eqn{\alpha_{c}}
#'   has the same prior distribution.
#'
#' @param psi Positive integer specifying the concentration parameter \eqn{\psi}
#'   of the Dirichlet prior distribution used for the cluster probabilities
#'   \eqn{\tau_{1}, \tau_{2}, \dots, \tau_{C}}, where \eqn{C} is the value of
#'   `n_clusters`. Defaults to `10L`. When `n_clusters = 1`, this argument is
#'   not used.
#'
#' @param kappa Hyperparameters of the truncated Beta prior used for error
#'   probability \eqn{\theta} in the Bernoulli error model. The prior has the
#'   form \eqn{\pi(\theta) = \theta^{\kappa_{1}} (1 - \theta)^{\kappa_{2}}}.
#'   Defaults to `c(1, 3)`, which means that the \eqn{\theta} is a priori
#'   expected to be closer to zero than to 0.5. See
#'   \insertCite{crispino2019}{BayesMallows} for details.
#'
#' @return An object of class `"BayesMallowsPriors"`, to be provided in the
#'   `priors` argument to [compute_mallows()], [compute_mallows_mixtures()], or
#'   [update_mallows()].
#' @export
#'
#' @references \insertAllCited{}
#'
#' @family preprocessing
#'
set_priors <- function(gamma = 1, lambda = 0.001, psi = 10, kappa = c(1, 3)) {
  stopifnot(length(kappa) == 2)
  validate_positive(gamma)
  validate_positive(lambda)
  validate_integer(psi)
  validate_positive(psi)
  validate_positive(kappa[[1]])
  validate_positive(kappa[[2]])

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsPriors"
  ret
}
