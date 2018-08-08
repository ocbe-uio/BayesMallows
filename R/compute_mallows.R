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
#' @param L Step size of the leap-and-shift proposal distribution. Defaults to
#'   NULL, which means that it is automatically set to n/5.
#' @param sd_alpha Standard deviation of the proposal distribution for alpha.
#' @param alpha_init Initial value of alpha.
#' @param alpha_jump How many times should we sample \code{rho} between each
#'   time we sample \code{alpha}. Setting \code{alpha_jump} to a high number can
#'   significantly speed up computation time, since we then do not have to do
#'   expensive computation of the partition function.
#' @param thinning Keep every \code{thinning} iteration of \code{rho}.
#'
#' @return A list of class BayesMallows.
#' @references \insertRef{vitelli2018}{BayesMallows}
#' @seealso \code{\link{assess_convergence}}, \code{\link{plot.BayesMallows}}.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' # Compute the posterior
#' model_fit <- compute_mallows(potato_weighing, "footrule", nmc = 1000)
#'
compute_mallows <- function(R,
                            metric = "footrule",
                            lambda = 0.1,
                            nmc = 3000,
                            L = NULL,
                            sd_alpha = 0.1,
                            alpha_init = 1,
                            alpha_jump = 1,
                            thinning = 1
                            ){

  # Check that all rows of R are proper permutations
  stopifnot(all(apply(R, 1, validate_permutation)))

  # Check that we do not jump over all alphas
  stopifnot(alpha_jump < nmc)

  # Check that we do not jump over all rhos
  stopifnot(thinning < nmc)

  # Find the number of items
  n_items <- ncol(R)

  # Set L if it is not alredy set.
  if(is.null(L)) L <- n_items / 5

  # Encoding missing values with -1
  R[is.na(R)] <- -1

  # Extract the right sequence of cardinalities, if relevant
  if(metric %in% c("footrule", "spearman")){
    # Extract the relevant rows from partition_function_data
    # Note that we need to evaluate the right-hand side, in particular metric,
    # to avoid confusion with columns of the tibble
    relevant_params <- dplyr::filter(partition_function_data,
                                     .data$num_items == !!n_items,
                                     .data$metric == !!metric
    )

    type <- dplyr::pull(relevant_params, type)

    if(type == "cardinalities") {
      cardinalities <- unlist(relevant_params$values)
      is_fit <- NULL
    } else if(type == "importance_sampling"){
      cardinalities <- NULL
      is_fit <- unlist(relevant_params$values)
    } else {
      stop("Precomputed partition function not available yet.")
    }

  } else if (metric %in% c("cayley", "hamming", "kendall")) {
    cardinalities <- NULL
    is_fit <- NULL
  } else {
    stop(paste("Unknown metric", metric))
  }

  # Transpose R to get samples along columns.
  fit <- run_mcmc(t(R), cardinalities = cardinalities, is_fit = is_fit,
                  metric = metric, lambda = lambda,
                  nmc = nmc, L = L, sd_alpha = sd_alpha,
                  alpha_init = alpha_init,
                  alpha_jump = alpha_jump, thinning = thinning)

  class(fit) <- "BayesMallows"

  return(fit)

}
