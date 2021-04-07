#' @title Metropolis-Hastings Alpha
#'  @decscription Function to perform Metropolis-Hastings for new alpha under the Mallows model
#'   @param alpha Numeric value og the scale parameter
#'   @param n_items Integer is the number of items in a ranking
#'   @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
#'   rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
#'   can be a vector.
#'   @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#'   @param rho Numeric vector specifying the consensus ranking
#'   @param logz_estimate Estimate of the partition function, computed with
#'   \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
#' @return \code{alpha} or \code{alpha_prime}: Numeric value to be used as the proposal of a new alpha
#'' @importFrom stats dexp rlnorm runif
metropolis_hastings_alpha <- function(alpha, n_items, rankings, metric, rho, logz_estimate){


  # select alpha_prime from log normal distribution with some variance lambda
  # log x ~ N(mu, sigma^2)
  # in C, this is written as
  # double alpha_proposal = std::exp(arma::randn<double>() * alpha_prop_sd + std::log(alpha_old));.
  exp_alpha_prime = rlnorm(1, meanlog = alpha, sdlog = 0.15) # 1
  alpha_prime = log(exp_alpha_prime)

  # evaluate the log-likelihood with current rankings
  mallows_loglik_prop = get_mallows_loglik(alpha = (alpha_prime - alpha), rho = rho, n_items = n_items, rankings = rankings,
                                           metric = metric)

  # evaluate the log estimate of the partition function for a particular value of alpha
  logz_alpha = get_partition_function(n_items = n_items, alpha = alpha, logz_estimate = logz_estimate, metric = metric)
  logz_alpha_prime = get_partition_function(n_items = n_items, alpha = alpha_prime, logz_estimate = logz_estimate, metric = metric)


  n_users = length(rankings)/n_items

  loga = n_users* (logz_alpha - logz_alpha_prime) + dexp(alpha_prime, log=TRUE) - dexp(alpha, log=TRUE) +
    alpha_prime - alpha + mallows_loglik_prop

  # determine whether to accept or reject proposed rho and return now consensus ranking
  p = runif(1, min = 0, max = 1)
  if(log(p) <= loga){
    return(alpha_prime)
  } else{
    return(alpha)
  }

}
