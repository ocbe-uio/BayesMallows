#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Alpha
//' @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!
//' @param alpha Numeric value og the scale parameter
//' @param n_items Integer is the number of items in a ranking
//' @param rankings the observed rankings, i.e, preference data
//' @details \code{rankings} is a matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
//' rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
//' can be a vector.
//' @param metric A character string specifying the distance metric to use in the
//' Bayesian Mallows Model. Available options are \code{"footrule"},
//' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//' \code{"ulam"}.
//' @param rho Numeric vector specifying the current consensus ranking
//' @param logz_estimate Estimate  grid of log of partition function, computed with
//' \code{\link{estimate_partition_function}} in the BayesMallow R package {estimate_partition_function}.
//' @return \code{alpha} or \code{alpha_prime}: Numeric value to be used as the proposal of a new alpha
//' @importFrom stats dexp rlnorm runif
//' @author Anja Stein
//' @example /inst/examples/metropolis_hastings_alpha.R
//'
//' @export
// [[Rcpp::export]]
// metropolis_hastings_alpha <- function(alpha, n_items, rankings, metric, rho, logz_estimate, alpha_prop_sd,
//                                              lambda, alpha_max){
//   alpha_prime <- exp(rnorm(n = 1, mean = 0, sd = 1) * alpha_prop_sd + log(alpha))

//   # Difference between current and proposed alpha
//   alpha_diff <- alpha - alpha_prime
//   mallows_loglik_prop <- get_mallows_loglik(
//     alpha = (alpha_prime - alpha), rho = rho, n = n_items, rankings = rankings,
//     metric = metric
//   )

//   # evaluate the log estimate of the partition function for a particular value of alpha
//   logz_alpha <- get_partition_function(
//     n_items = n_items, alpha = alpha,
//     logz_estimate = logz_estimate, metric = metric
//   )
//   logz_alpha_prime <- get_partition_function(
//     n_items = n_items, alpha = alpha_prime,
//     logz_estimate = logz_estimate, metric = metric
//   )

//   obs_freq <- length(rankings) / n_items

//   # Compute the Metropolis-Hastings ratio
//   loga <- mallows_loglik_prop +
//     lambda * alpha_diff +
//     obs_freq * (logz_alpha - logz_alpha_prime) +
//     log(alpha_prime) - log(alpha)

//   # determine whether to accept or reject proposed rho and return now consensus
//   # ranking
//   p <- runif(1, min = 0, max = 1)
//   if(log(p) <= loga && alpha_prime < alpha_max) {
//     return(alpha_prime)
//   } else {
//     return(alpha)
//   }
// }