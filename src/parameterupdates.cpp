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
