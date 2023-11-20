#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"
#include "partitionfuns.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Initialize latent ranks as provided by rho_init, or randomly:
arma::mat initialize_rho(int n_items, int n_cols,
                         Rcpp::Nullable<arma::mat> rho_init){
  if(rho_init.isNotNull()){
    return repmat(Rcpp::as<mat>(rho_init), 1, n_cols);
  } else {
    mat rho_samples(n_items, n_cols);
    for (int i = 0; i < n_cols; ++i) {
      rho_samples.col(i) = randperm<vec>(n_items) + 1;
    }
    return rho_samples;
  }
}

double update_alpha(
                  const double& alpha_old,
                  const mat& rankings,
                  const vec& observation_frequency,
                  const vec& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::List& logz_list) {
  // Set the number of assessors. Not using the variable from run_mcmc because
  // here we want the number of assessors in this cluster
  //int n_assessors = rankings.n_cols;
  int n_items = rho_old.n_elem;

  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                              std::log(alpha_old));

  double rank_dist = rank_dist_sum(rankings, rho_old, metric, observation_frequency);

  // Difference between current and proposed alpha
  double alpha_diff = alpha_old - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    sum(observation_frequency) * (
        get_partition_function(n_items, alpha_old, logz_list, metric) -
          get_partition_function(n_items, alpha_proposal, logz_list, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  if(ratio > u){
    return alpha_proposal;
  } else {
    return alpha_old;
  }
}


vec make_new_rho(vec current_rho, const mat& rankings, double alpha_old, int leap_size, std::string metric,
                 vec observation_frequency) {

  // Sample a rank proposal
  vec rho_proposal;
  uvec indices;
  double prob_backward, prob_forward;
  int n_items = current_rho.n_elem;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 current_rho, leap_size, !((metric == "cayley") || (metric == "ulam")));

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric, observation_frequency);
  double dist_old = rank_dist_sum(rankings.rows(indices), current_rho(indices), metric, observation_frequency);

  // Metropolis-Hastings ratio
  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  if(ratio > u){
    return(rho_proposal);
  } else {
    return(current_rho);
  }

}
