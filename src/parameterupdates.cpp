#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

void update_alpha(arma::mat& alpha,
                  arma::mat& alpha_acceptance,
                  arma::vec& alpha_old,
                  const arma::mat& R,
                  int& alpha_index,
                  int& cluster_index,
                  const arma::mat& rho_old,
                  const double& sd_alpha,
                  const std::string& metric,
                  const double& lambda,
                  const int& n_items,
                  Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  Rcpp::Nullable<arma::vec> is_fit = R_NilValue) {


  // Set the number of assessors. Not using the variable from run_mcmc because
  // here we want the number of assessors in this cluster #cluster_index
  int n_assessors = R.n_cols;

  // Sample an alpha proposal
  double alpha_proposal = exp(arma::randn<double>() * sd_alpha +
  log(alpha(alpha_index - 1, cluster_index)));

  double rank_dist = rank_dist_matrix(R, rho_old.col(cluster_index), metric);

  // Difference between current and proposed alpha
  double alpha_diff = alpha(alpha_index - 1, cluster_index) - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    n_assessors * (
        get_partition_function(n_items, alpha(alpha_index - 1, cluster_index),
                               cardinalities, is_fit, metric) -
          get_partition_function(n_items, alpha_proposal, cardinalities, is_fit, metric)
    ) + log(alpha_proposal) - log(alpha(alpha_index - 1, cluster_index));

  // Draw a uniform random number
  double u = log(arma::randu<double>());




  if(ratio > u){
    alpha(alpha_index, cluster_index) = alpha_proposal;
    alpha_acceptance(alpha_index, cluster_index) = 1;
  } else {
    alpha(alpha_index, cluster_index) = alpha(alpha_index - 1, cluster_index);
    alpha_acceptance(alpha_index, cluster_index) = 0;
  }


  alpha_old(cluster_index) = alpha(alpha_index, cluster_index);

}


void update_rho(arma::cube& rho, arma::mat& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& thinning,
                const double& alpha_old, const int& L, const arma::mat& R,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices) {


  // Sample a rank proposal
  Rcpp::List ls_proposal = leap_and_shift(rho_old.col(cluster_index), L);

  // Save some of the variables
  arma::vec rho_proposal = ls_proposal["proposal"];
  arma::uvec indices = ls_proposal["indices"];
  double prob_backward = ls_proposal["prob_backward"];
  double prob_forward = ls_proposal["prob_forward"];

  // Compute the distances to current and proposed ranks
  double dist_new = rank_dist_matrix(R.rows(indices), rho_proposal(indices), metric);

  double dist_old = rank_dist_matrix(R.rows(indices),
                                     rho_old.submat(indices, arma::regspace<arma::uvec>(cluster_index, cluster_index)),
                                     metric);


  // Metropolis-Hastings ratio
  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    log(prob_backward) - log(prob_forward);

  // Draw a uniform random number
  double u = log(arma::randu<double>());

  if(ratio > u){
    rho_old.col(cluster_index) = rho_proposal;
    rho_acceptance.col(cluster_index).row(t) = 1;
  } else {
    rho_acceptance.col(cluster_index).row(t) = 0;
  }

  // Save rho if appropriate
  if(t % thinning == 0){
    if(cluster_index == 0) ++rho_index;
    rho.slice(cluster_index).col(rho_index) = rho_old.col(cluster_index);
  }

}


void update_cluster_labels(
    arma::umat& cluster_indicator,
    const arma::mat& rho_old,
    const arma::mat& R,
    const arma::mat& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const int& n_assessors,
    const int& n_clusters,
    const int& t,
    const std::string& metric,
    Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    Rcpp::Nullable<arma::vec> is_fit = R_NilValue
){

  // Matrix to hold assignment probabilities
  arma::mat assignment_prob(n_clusters, n_assessors);

  for(int assessor_index = 0; assessor_index < n_assessors; ++assessor_index){
    for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
      assignment_prob(cluster_index, assessor_index) =
        cluster_probs(cluster_index, t) *
          exp(-alpha_old(cluster_index) / n_items *
          get_rank_distance(R.col(assessor_index), rho_old.col(cluster_index), metric)) /
            exp(get_partition_function(n_items, alpha_old(cluster_index), cardinalities, is_fit, metric));
    }

  }

  // Normalise the assignment probabilities, to unit L1 norm for each assessor (column)
  assignment_prob = arma::normalise(assignment_prob, 1, 0);

  for(int assessor_index = 0; assessor_index < n_assessors; ++assessor_index){
    // Draw a uniform random number
    double u = arma::randu();

    // Pick out the column vector of probabilities for this assessor
    // This corresponds to aa in JMLR code
    arma::vec assessor_prob_vector = assignment_prob.col(assessor_index);

    // Define helper variables
    bool cluster_assigned = false;
    int cluster_index = 0;
    double sum = 0;

    while(!cluster_assigned){
      if(u < (arma::as_scalar(assessor_prob_vector(cluster_index)) / (1 - sum))) {
        cluster_assigned = true;
      } else {
        sum += arma::as_scalar(assessor_prob_vector(cluster_index));
        ++cluster_index;
      }
    }

    cluster_indicator(t, assessor_index) = cluster_index;
  }

}
