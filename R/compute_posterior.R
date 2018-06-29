#' Compute the posterior distribution for the Mallows model.
#'
#' @param R A matrix of ranks, N x n. This matrix is transposed before calling
#'   the underlying C++ function.
#' @param metric The metric to use.
#'
#' @return TBD
#' @export
#'
#' @examples
#' model_fit <- compute_posterior(potato_weighing, "footrule")
compute_posterior <- function(R, metric = "footrule", lambda = 0.1,
                              nmc = 1000, L = 1, sd_alpha = 0.1,
                              alpha_init = 5){
  # Find the number of items
  n <- ncol(R)

  # Extract the right sequence of cardinalities, if relevant
  if(metric == "footrule"){
    cardinalities <- footrule_sequence[[n]]
  } else if (metric == "spearman") {
    cardinalities <- spearman_sequence[[n]]
  } else {
    cardinalities <- 0
  }

  # Transpose R to get samples along columns.
  run_mcmc(t(R), cardinalities, metric = metric, lambda = lambda,
           nmc = nmc, L = L, sd_alpha = sd_alpha, alpha_init = alpha_init)

}
