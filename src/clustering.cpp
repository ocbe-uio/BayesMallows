#include "RcppArmadillo.h"
#include <cmath>
#include "partitionfuns.h"
#include "misc.h"

// [[Rcpp::depends(RcppArmadillo)]]

void update_cluster_labels(
    arma::uvec& current_cluster_assignment,
    const arma::mat& dist_mat,
    const arma::mat& rho_old,
    const arma::mat& rankings,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_old,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue
){
  int n_assessors = current_cluster_assignment.n_elem;
  int n_items = rho_old.n_rows;
  int n_clusters = rho_old.n_cols;

  arma::mat assignment_prob(n_assessors, n_clusters);
  for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(cluster_index) = std::log(cluster_probs(cluster_index)) -
      alpha_old(cluster_index) / n_items * dist_mat.col(cluster_index) -
      get_partition_function(n_items, alpha_old(cluster_index), cardinalities, logz_estimate, metric);
  }

  for(int assessor_index = 0; assessor_index < n_assessors; ++assessor_index){
    // Exponentiate to get unnormalized prob relative to max
    arma::rowvec probs = arma::exp(assignment_prob.row(assessor_index) -
      arma::max(assignment_prob.row(assessor_index)));

    // Normalize with 1-norm
    probs = arma::normalise(probs, 1);

    assignment_prob.row(assessor_index) = probs;
    int cluster = sample_int(assignment_prob.row(assessor_index));

    // Assign the cluster indicator
    current_cluster_assignment(assessor_index) = cluster;

  }
}

void update_cluster_probs(
    arma::mat& cluster_probs,
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi,
    const int& t
){
  // Update cluster probabilities, conjugate model
  arma::uvec tau_k(n_clusters);
  for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
    // Find the parameter for this cluster
    tau_k(cluster_index) = arma::sum(current_cluster_assignment == cluster_index) + psi;

    // If there are no assessors in the cluster,
    // Save a sample from the gamma distribution
    cluster_probs(cluster_index, t) = arma::randg<double>(arma::distr_param(tau_k(cluster_index), 1.0));
  }

  // Finally, normalize cluster_probs. Must specify that it should have unit 1-norm,
  // 2-norm is default!!
  // cluster_probs.col(t) now comes from Dirichlet(tau_k(0), ..., tau_k(n_clusters))
  cluster_probs.col(t) = arma::normalise(cluster_probs.col(t), 1);
}


void update_wcd(arma::mat& within_cluster_distance, const arma::uvec& current_cluster_assignment,
                const arma::mat& dist_mat, const int& n_clusters, const int& t){

  arma::vec dist_vec;

  for(int i = 0; i < n_clusters; ++i){
    // Find the assessors in cluster i
    arma::uvec assessor_inds = arma::find(current_cluster_assignment == i);

    // Extract the distances from R to rho_i
    dist_vec = dist_mat.col(i);

    // Extract the assessors in cluster i
    dist_vec = dist_vec(assessor_inds);

    // Sum and save
    within_cluster_distance(i, t) = arma::sum(dist_vec);
  }
}
