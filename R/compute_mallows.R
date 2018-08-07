#' Compute the posterior distribution for the Mallows model.
#'
#' @param R A matrix of ranked items. See \code{\link{create_ranking}} if you
#'   have an ordered set of items that needs to be converted to rankings.
#' @param metric The distance metric to use. Currently available ones are
#'   \code{"footrule"}, \code{"spearman"}, \code{"kendall"}, \code{"cayley"},
#'   and \code{"hamming"}.
#' @param lambda Parameter for the prior distribution of \code{alpha}. Defaults
#'   to 0.1.
#' @param nmc Number of Monte Carlo samples to keep.
#' @param burnin Number of samples to discard as burn-in before collecting
#'   \code{nmc} samples. See \code{\link{assess_convergence}} for a function
#'   that can be used to pick the right number of burn-in samples to discard.
#' @param L Step size of the leap-and-shift proposal distribution. Defaults to
#'   NULL, which means that it is automatically set to n/5.
#' @param sd_alpha Standard deviation of the proposal distribution for alpha.
#' @param alpha_init Initial value of alpha.
#' @param alpha_jump How many times should we sample \code{rho} between each
#'   time we sample \code{alpha}. Setting \code{alpha_jump} to a high number can
#'   significantly speed up computation time, since we then do not have to do
#'   expensive computation of the partition function.
#'
#' @return A list of class BayesMallows.
#' @details  It is usually a good idea to first use
#'   \code{\link{assess_convergence}} to determine the algorithm parameters.
#' @references \insertRef{vitelli2018}{BayesMallows}
#' @seealso \code{\link{assess_convergence}}, \code{\link{plot.BayesMallows}}.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # Compute the posterior
#' model_fit <- compute_mallows(potato_weighing, "footrule", nmc = 10000, burnin = 5000)
#' # Plot the posterior histogram for alpha
#' plot(model_fit, type = "alpha", bins = 50)
#' # Plot the posterior histogram for potatoes 1, 2, and 3.
#' plot(model_fit, type = "rho", items = 1:3)
#'
compute_mallows <- function(R, metric = "footrule", lambda = 0.1,
                              nmc = 3000, burnin = 2000, L = NULL, sd_alpha = 0.1,
                              alpha_init = 0.001, alpha_jump = 1){

  # Check that we have more samples than we throw away
  stopifnot(nmc > burnin)

  # Check that all rows of R are proper permutations
  stopifnot(all(apply(R, 1, validate_permutation)))

  # Check that we do not jump over all alphas
  stopifnot(alpha_jump < nmc)

  # Find the number of items
  n_items <- ncol(R)

  # Set L if it is not alredy set.
  if(is.null(L)) L <- n_items / 5

  ## Temporary!
  is_fit <- NULL

  # Extract the right sequence of cardinalities, if relevant
  if(metric %in% c("footrule", "spearman")){
    # Extract the relevant rows from partition_function_data
    # Note that we need to evaluate the right-hand side, in particular metric,
    # to avoid confusion with columns of the tibble
    relevant_params <- dplyr::filter(partition_function_data,
                                     .data$num_items == !!n_items,
                                     .data$metric == !!metric,
                                     .data$type == "cardinalities" # TEMPORARY!!
    )

    cardinalities <- get_cardinalities(relevant_params)
  } else if (metric %in% c("cayley", "hamming", "kendall")) {
    cardinalities <- NULL
  } else {
    stop(paste("Unknown metric", metric))
  }

  # Transpose R to get samples along columns.
  fit <- run_mcmc(t(R), cardinalities = cardinalities, is_fit = is_fit, metric = metric, lambda = lambda,
           nmc = nmc, L = L, sd_alpha = sd_alpha, alpha_init = alpha_init,
           alpha_jump = alpha_jump)

  # Finally subset so that we only get values after burnin
  fit$rho <- fit$rho[, (burnin + 1) : nmc]
  fit$rho_acceptance <- fit$rho_acceptance[(burnin + 1) : nmc]

  fit$alpha <- fit$alpha[(floor((burnin + 1) / alpha_jump)) : (floor(nmc / alpha_jump))]
  fit$alpha_acceptance <- fit$alpha_acceptance[(floor((burnin + 1) / alpha_jump)) : (floor(nmc / alpha_jump))]

  class(fit) <- "BayesMallows"

  return(fit)

}
