#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


void update_cluster_labels(
    arma::umat& cluster_assignment,
    const arma::mat& dist_mat,
    const arma::mat& rho_old,
    const arma::mat& rankings,
    const arma::mat& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const int& n_assessors,
    const int& n_clusters,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> is_fit = R_NilValue
){

  arma::mat assignment_prob(n_assessors, n_clusters);

  for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
    assignment_prob.col(cluster_index) = cluster_probs(cluster_index, t) *
      arma::exp(-alpha_old(cluster_index) / n_items * dist_mat.col(cluster_index)) /
        exp(get_partition_function(n_items, alpha_old(cluster_index), cardinalities, is_fit, metric));
  }

  // Normalize the probabilities, to unit 1-norm per row (assessor)
  assignment_prob = arma::normalise(assignment_prob, 1, 1);

  for(int assessor_index = 0; assessor_index < n_assessors; ++assessor_index){

    int cluster = sample_int(assignment_prob.row(assessor_index));

    // Assign the cluster indicator
    cluster_assignment(assessor_index, t) = cluster;


  }

}

void update_cluster_probs(
    arma::mat& cluster_probs,
    const arma::umat& cluster_assignment,
    const int& n_clusters,
    const int& psi,
    const int& t
){
  // Update cluster probabilities, conjugate model
  arma::uvec tau_k(n_clusters);
  for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
    // Find the parameter for this cluster
    tau_k(cluster_index) = arma::sum(cluster_assignment.col(t - 1) == cluster_index) + psi;

    // If there are no assessors in the cluster,
    // Save the a draw from the gamma distribution
    cluster_probs(cluster_index, t) = arma::randg<double>(arma::distr_param(tau_k(cluster_index), 1.0));

  }

  // Finally, normalize cluster_probs. Must specify that it should have unit 1-norm,
  // 2-norm is default!!
  // cluster_probs.col(t) now comes from Dirichlet(tau_k(0), ..., tau_k(n_clusters))
  cluster_probs.col(t) = arma::normalise(cluster_probs.col(t), 1);
}


void update_wcd(arma::mat& within_cluster_distance, const arma::uvec& cluster_assignment,
                const arma::mat& dist_mat, const int& n_clusters, const int& t){

  arma::vec dist_vec;

  for(int i = 0; i < n_clusters; ++i){
    // Find the assessors in cluster i
    arma::uvec assessor_inds = arma::find(cluster_assignment == i);

    // Extract the distances from R to rho_i
    dist_vec = dist_mat.col(i);

    // Extract the assessors in cluster i
    dist_vec = dist_vec(assessor_inds);

    // Sum and save
    within_cluster_distance(i, t) = arma::sum(dist_vec);
  }
}
