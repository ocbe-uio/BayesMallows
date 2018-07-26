#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

void update_alpha(arma::vec& alpha, arma::vec& alpha_acceptance,
                  double& alpha_old,
                  const arma::mat& R, int& alpha_index,
                  const arma::vec& rho_old,
                  const double& sd_alpha, const std::string& metric,
                  const double& lambda, const int& n, const int& N,
                  Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  Rcpp::Nullable<arma::vec> is_fit = R_NilValue) {

  // Increment to the index we are going to update
  ++alpha_index;

  // Sample an alpha proposal (normal on the log scale)
  double alpha_proposal = exp(arma::randn<double>() * sd_alpha +
  log(alpha(alpha_index - 1)));

  double rank_dist = rank_dist_matrix(R, rho_old, metric);

  // Difference between current and proposed alpha
  double alpha_diff = alpha(alpha_index - 1) - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio = (alpha(alpha_index - 1) - alpha_proposal) / n * rank_dist +
    lambda * alpha_diff +
    N * (
        get_partition_function(n, alpha(alpha_index - 1), cardinalities, is_fit, metric) -
          get_partition_function(n, alpha_proposal, cardinalities, is_fit, metric)
    ) + log(alpha_proposal) - log(alpha(alpha_index - 1));

  // Draw a uniform random number
  double u = log(arma::randu<double>());


  if(ratio > u){
    alpha(alpha_index) = alpha_proposal;
    alpha_acceptance(alpha_index) = 1;
  } else {
    alpha(alpha_index) = alpha(alpha_index - 1);
    alpha_acceptance(alpha_index) = 0;
  }

  // Next time around, this is alpha_old
  // The sampling of rho below should also
  // use this alpha_old, which then in the first
  // step is actually the newest alpha
  alpha_old = alpha(alpha_index);
}


void update_rho(arma::mat& rho, arma::vec& rho_acceptance, arma::vec& rho_old,
                const double& alpha_old, const int& L, const arma::mat& R,
                const std::string& metric, const int& n, const int& t) {
  // Sample a rank proposal
  Rcpp::List ls_proposal = leap_and_shift(rho_old, L);

  // Save some of the variables
  arma::vec rho_proposal = ls_proposal["proposal"];
  arma::uvec indices = ls_proposal["indices"];
  double prob_backward = ls_proposal["prob_backward"];
  double prob_forward = ls_proposal["prob_forward"];

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_matrix(R.rows(indices), rho_proposal(indices), metric);
  double dist_old = rank_dist_matrix(R.rows(indices), rho_old(indices), metric);

  // Metropolis-Hastings ratio
  double ratio = - alpha_old / n * (dist_new - dist_old) +
    log(prob_backward) - log(prob_forward);

  // Draw a uniform random number
  double u = log(arma::randu<double>());

  if(ratio > u){
    rho.col(t) = rho_proposal;
    rho_acceptance(t) = 1;
  } else {
    rho.col(t) = rho.col(t - 1);
    rho_acceptance(t) = 0;
  }

  rho_old = rho.col(t - 1);
}
