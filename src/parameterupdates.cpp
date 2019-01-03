#include "RcppArmadillo.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"


// [[Rcpp::depends(RcppArmadillo)]]


double rtruncbeta(int shape1, int shape2, double trunc = 1) {
  int i = 0;
  double x;
  while(i < 1000){
    x = arma::chi2rnd(2 * shape1);
    x = x / (x + arma::chi2rnd(2 * shape2));

    if(x < trunc) break;
    ++i;
  }
  return x;
}


void update_alpha(arma::mat& alpha,
                  arma::vec& alpha_acceptance,
                  arma::vec& alpha_old,
                  const arma::mat& rankings,
                  int& alpha_index,
                  int& cluster_index,
                  const arma::mat& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const int& n_items,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue) {


  // Set the number of assessors. Not using the variable from run_mcmc because
  // here we want the number of assessors in this cluster #cluster_index
  int n_assessors = rankings.n_cols;

  // Sample an alpha proposal
  double alpha_proposal = std::exp(arma::randn<double>() * alpha_prop_sd +
                              std::log(alpha_old(cluster_index)));

  double rank_dist = rank_dist_sum(rankings, rho_old.col(cluster_index), metric);

  // Difference between current and proposed alpha
  double alpha_diff = alpha_old(cluster_index) - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    n_assessors * (
        get_partition_function(n_items, alpha_old(cluster_index), cardinalities, logz_estimate, metric) -
          get_partition_function(n_items, alpha_proposal, cardinalities, logz_estimate, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old(cluster_index));

  // Draw a uniform random number
  double u = std::log(arma::randu<double>());

  if(ratio > u){
    alpha(cluster_index, alpha_index) = alpha_proposal;
    ++alpha_acceptance(cluster_index);
  } else {
    alpha(cluster_index, alpha_index) = alpha_old(cluster_index);
  }


  alpha_old(cluster_index) = alpha(cluster_index, alpha_index);

}


void update_rho(arma::cube& rho, arma::vec& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& rho_thinning,
                const double& alpha_old, const int& leap_size, const arma::mat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices, bool& rho_accepted) {

  arma::vec rho_cluster = rho_old.col(cluster_index);

  // Sample a rank proposal
  arma::vec rho_proposal;
  arma::uvec indices;
  double prob_backward, prob_forward;

  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho_cluster, leap_size);

  // These distances do not work with the computational shortcut
  if((metric == "cayley") | (metric == "ulam")){
    indices = arma::regspace<arma::uvec>(0, n_items - 1);
  }

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric);
  double dist_old = rank_dist_sum(rankings.rows(indices), rho_cluster(indices), metric);


  // Metropolis-Hastings ratio
  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  // Draw a uniform random number
  double u = std::log(arma::randu<double>());

  if(ratio > u){
    rho_old.col(cluster_index) = rho_proposal;
    ++rho_acceptance(cluster_index);
    rho_accepted = true;
  } else {
    rho_accepted = false;
  }

  // Save rho if appropriate
  if(t % rho_thinning == 0){
    if(cluster_index == 0) ++rho_index;
    rho.slice(rho_index).col(cluster_index) = rho_old.col(cluster_index);
  }

}



