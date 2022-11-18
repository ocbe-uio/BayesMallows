#' @title Plot SMC Posterior Distributions
#' @description Plot posterior distributions of SMC-Mallow parameters.
#' @param x An object of type \code{SMC-Mallows}, returned for example from
#' [smc_mallows_new_users()].
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See
#' [assess_convergence()]
#' @param parameter Character string defining the parameter to plot. Available
#' options are \code{"alpha"} and \code{"rho"}.
#' @param time Integer determining the update slice to plot
#' @param C Number of cluster
#' @param colnames A vector of item names. If NULL, generic names are generated
#' for the items in the ranking.
#' @param items Either a vector of item names, or a vector of indices. If NULL,
#' five items are selected randomly.
#' @param ... Other arguments passed to [base::plot()] (not used).
#' @return A plot of the posterior distributions
#' @author Waldir Leoncio
#' @export
#' @example /inst/examples/plot.SMCMallows_example.R
plot.SMCMallows <- function(x, nmc = nrow(x$rho_samples[, 1, ]), burnin = 0,
  parameter = "alpha", time = ncol(x$rho_samples[, 1, ]), C = 1,
  colnames = NULL, items = NULL, ...) {

  if (parameter == "alpha") {
    output <- x$alpha_samples[, time]
    plot_alpha_posterior(output, nmc, burnin)
  } else if (parameter == "rho") {
    output <- x$rho_samples[, , time]
    plot_rho_posterior(output, nmc, burnin, C, colnames, items)
  } else {
    stop("parameter must be either 'alpha' or 'rho'.")
  }
}
