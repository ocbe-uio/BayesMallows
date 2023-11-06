#include <RcppArmadillo.h>
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

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
