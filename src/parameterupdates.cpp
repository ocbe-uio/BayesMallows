#include "RcppArmadillo.h"
#include "leapandshift.h"
#include "distances.h"
#include "partitionfuns.h"
#include <cmath>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Initialize latent ranks as provided by rho_init, or randomly:
mat initialize_rho(Rcpp::Nullable<mat> rho_init, int n_items, int n_clusters){
  if(rho_init.isNotNull()){
    return arma::repmat(Rcpp::as<mat>(rho_init), 1, n_clusters);
  } else {
    return arma::shuffle(arma::repmat(regspace<mat>(1, 1, n_items), 1, n_clusters));
  }
}

double update_alpha(vec& alpha_acceptance,
                  const double& alpha_old,
                  const mat& rankings,
                  const vec& obs_freq,
                  const int& cluster_index,
                  const vec& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::Nullable<vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<vec> logz_estimate = R_NilValue,
                  double alpha_max = 1e6) {


  // Set the number of assessors. Not using the variable from run_mcmc because
  // here we want the number of assessors in this cluster
  //int n_assessors = rankings.n_cols;
  int n_items = rho_old.n_elem;

  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                              std::log(alpha_old));

  double rank_dist = rank_dist_sum(rankings, rho_old, metric, obs_freq);


  // Difference between current and proposed alpha
  double alpha_diff = alpha_old - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    sum(obs_freq) * (
        get_partition_function(n_items, alpha_old, cardinalities, logz_estimate, metric) -
          get_partition_function(n_items, alpha_proposal, cardinalities, logz_estimate, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  if(ratio > u && alpha_proposal < alpha_max){
    ++alpha_acceptance(cluster_index);
    return alpha_proposal;
  } else {
    return alpha_old;
  }
}


void update_rho(cube& rho, vec& rho_acceptance, mat& rho_old,
                int& rho_index, const int& cluster_index, const int& rho_thinning,
                const double& alpha_old, const int& leap_size, const mat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const uvec& element_indices, const vec& obs_freq) {

  vec rho_cluster = rho_old.col(cluster_index);

  // Sample a rank proposal
  vec rho_proposal;
  uvec indices;
  double prob_backward, prob_forward;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho_cluster, leap_size, !((metric == "cayley") || (metric == "ulam")));


  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric, obs_freq);
  double dist_old = rank_dist_sum(rankings.rows(indices), rho_cluster(indices), metric, obs_freq);

  // Metropolis-Hastings ratio
  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  if(ratio > u){
    rho_old.col(cluster_index) = rho_proposal;
    ++rho_acceptance(cluster_index);
  }

  // Save rho if appropriate
  if(t % rho_thinning == 0){
    if(cluster_index == 0) ++rho_index;
    rho.slice(rho_index).col(cluster_index) = rho_old.col(cluster_index);
  }

}



