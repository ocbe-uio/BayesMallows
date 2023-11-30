#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"
#include "partitionfuns.h"
#include "parameterupdates.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


AlphaRatio make_new_alpha(const double& alpha_old, const vec& rho_old,
                          const double& alpha_prop_sd, const std::string& metric,
                          const Rcpp::List& logz_list,
                          const arma::mat& rankings,
                          const arma::vec& observation_frequency,
                          const double& n_items,
                          const Priors& priors) {
  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                                   std::log(alpha_old));

  double rank_dist = rank_dist_sum(
    rankings, rho_old, metric, observation_frequency);

  // Difference between current and proposed alpha
  double alpha_diff = alpha_old - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    priors.lambda * alpha_diff +
    sum(observation_frequency) * (
        get_partition_function(n_items, alpha_old, logz_list, metric) -
          get_partition_function(n_items, alpha_proposal, logz_list, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  return AlphaRatio{alpha_proposal, ratio > u};
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
