#include <RcppArmadillo.h>
#include "smc.h"
#include "partitionfuns.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double metropolis_hastings_alpha(
  const double alpha_old,
  const int n_items,
  const arma::mat rankings,
  const arma::vec rho,
  const Rcpp::Nullable<arma::vec> logz_estimate,
  const Rcpp::Nullable<arma::vec> cardinalities,
  const std::string metric = "footrule",
  const double alpha_prop_sd = 0.5,
  const double lambda = 0.1
) {

  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                                   std::log(alpha_old));

  int n_assessors = rankings.n_cols;
  vec obs_freq = ones(n_assessors);


  double rank_dist = rank_dist_sum(rankings, rho, metric, obs_freq);

  double alpha_diff = alpha_old - alpha_proposal;

  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    sum(obs_freq) * (
        get_partition_function(n_items, alpha_old, cardinalities, logz_estimate, metric) -
          get_partition_function(n_items, alpha_proposal, cardinalities, logz_estimate, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old);

  double u = std::log(randu<double>());

  if(ratio > u){
    return alpha_proposal;
  } else {
    return alpha_old;
  }

}
