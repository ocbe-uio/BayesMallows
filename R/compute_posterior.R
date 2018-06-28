#' Compute the posterior distribution for the Mallows model.
#'
#' @param R A matrix of ranks.
#' @param metric The metric to use.
#'
#' @return TBD
#' @export
#'
#' @examples
#' compute_posterior(potato_weighing, "footrule")
compute_posterior <- function(R, metric){
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

  run_mcmc(as.numeric(R, cardinalities, metric)

}
