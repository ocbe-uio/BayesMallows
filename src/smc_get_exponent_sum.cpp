#include <RcppArmadillo.h>
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Get exponent in Mallows log-likelihood
//' @description Calculates the exponent Mallows log-likelihood given a set of rankings
//' and a given rank sequence.
//' @param alpha Numeric value of the scale parameter
//' @param rho A ranking sequence
//' @param n_items Integer is the number of items in a ranking
//' A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
//' rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
//' can be a vector.
//' @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
//' rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
//' can be a vector.
//' @param metric Character string specifying the distance measure to use.
//' Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"},
//' \code{"ulam"}, \code{"footrule"} and \code{"spearman"}.
//' @return Exponent in the Mallows log likelihood. Note that it does not include
//' the partition function, and since the partition function depends on \code{alpha},
//' this is not a likelihood per se.
//' @noRd
//' @examples
//' set.seed(101)
//' rho <- t(c(1, 2, 3, 4, 5, 6))
//' alpha <- 2
//' metric <- "footrule"
//' n_items <- 6
//' get_exponent_sum(
//'   alpha = alpha, rho = rho, n_items = length(rho), rankings = rho,
//'   metric = metric
//' )
//'
//' # return 0 because you are comparing the consensus ranking with itself
//' # if you change alpha or metric, then the result shall remain as 0
//'
//' rankings <- sample_mallows(
//'   rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
//' )
//'
//' # depending on your seed, you will get a different collection of rankings in R and C++
//'
//' get_exponent_sum(
//'   alpha = alpha, rho = rho,  n_items = n_items, rankings = rankings ,
//'   metric = metric
//' )
// [[Rcpp::export]]
double get_exponent_sum(
  const double alpha,
  const arma::vec rho,
  const int n_items,
  arma::mat rankings,
  const std::string metric = "footrule"
) {

  vec obs_freq = ones(rankings.n_cols);
  double sum_distance = rank_dist_sum(rankings, rho, metric, obs_freq);
  double mallows_loglik = -alpha / n_items * sum_distance;
  return(mallows_loglik);
}
