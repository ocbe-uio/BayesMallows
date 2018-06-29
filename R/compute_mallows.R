#' Compute the posterior distribution for the Mallows model.
#'
#' @param R A matrix of ranked items.
#' @param metric The distance metric to use.
#' @param lambda Parameter for the prior distribution of \code{alpha}. Defaults to 0.1.
#' @param nmc Number of Monte Carlo samples to keep.
#' @param burnin Number of samples to discard as burn-in before collecting \code{nmc} samples.
#' @param L Step size of the leap-and-shift proposal distribution.
#' @param sd_alpha Standard deviation of the proposal distribution for alpha.
#' @param alpha_init Initial value of alpha.
#'
#' @return TBD
#' @details  It is usually a
#' good idea to first use \code{\link{assess_convergence}} to determine the
#' algorithm parameters.
#' @seealso \code{\link{assess_convergence}}, \code{\link{plot.BayesMallows}}.
#' @export
#'
#' @examples
#' # Compute the posterior
#' model_fit <- compute_mallows(potato_weighing, "footrule", nmc = 10000, burnin = 5000)
#' # Plot the posterior histogram for alpha
#' plot(model_fit, type = "alpha", bins = 50)
compute_mallows <- function(R, metric = "footrule", lambda = 0.1,
                              nmc = 3000, burnin = 2000, L = 1, sd_alpha = 0.1,
                              alpha_init = 5){

  # Check that we have more samples than we throw away
  stopifnot(nmc > burnin)

  # Find the number of items
  n <- ncol(R)

  # Extract the right sequence of cardinalities, if relevant
  if(metric == "footrule"){
    if(n > length(footrule_sequence)) {
      stop("At the moment, at most", length(footrule_sequence), "items can be analyzed.")
    }
    cardinalities <- footrule_sequence[[n]]
  } else if (metric == "spearman") {
    stop("Spearman distance is not implemented yet.")
    cardinalities <- spearman_sequence[[n]]
  } else {
    stop("Only footrule is implemented at the moment.")
    cardinalities <- 0
  }

  # Transpose R to get samples along columns.
  fit <- run_mcmc(t(R), cardinalities, metric = metric, lambda = lambda,
           nmc = nmc, L = L, sd_alpha = sd_alpha, alpha_init = alpha_init)

  # Finally subset so that we only get values after burnin
  s <- c("rho_acceptance", "alpha", "alpha_acceptance")
  fit[s] <- lapply(fit[s], `[`, seq(from = burnin + 1, to = nmc, by = 1))
  fit$rho <- fit$rho[, (burnin + 1) : nmc]
  class(fit) <- "BayesMallows"

  return(fit)

}
