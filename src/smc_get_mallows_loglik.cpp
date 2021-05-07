#include "RcppArmadillo.h"
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Get Mallows log-likelihood (CPP version)
//' @description Calculates the Mallows log-likelihood given a set of rankings and a given rank sequence
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
//' @return Mallows log-likelihood
//' @export
//' @author Anja Stein
//' @examples
//' set.seed(101)
//' rho <- c(1,2,3,4,5,6)
//' alpha <- 2
//' metric <- "footrule"
//' n_items <- 6
//' get_mallows_loglik(
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
//' get_mallows_loglik(
//'   alpha = alpha, rho = rho,  n_items = n_items, rankings = rankings ,
//'   metric = metric
//' )
// [[Rcpp::export]]
double get_mallows_loglik_CPP(
  double alpha,
  arma::vec rho,
  int n_items,
  arma::mat rankings,
  std::string metric
) {
  double sum_distance = 0;
  arma::uword num_rankings = rankings.n_rows;

  /* calculate the sum of the distances ------------------- */
  if (num_rankings == 1) {
    sum_distance = sum_distance + get_rank_distance(rho, rankings, metric);
  } else {
    for (int jj = 0; jj < num_rankings; ++jj) {
      sum_distance = sum_distance + get_rank_distance(rho, rankings.row(jj), metric);
    }
  }
  double mallows_loglik = -alpha / n_items * sum_distance;
  return(mallows_loglik);
}
