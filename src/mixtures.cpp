#include "RcppArmadillo.h"
#include <cmath>
#include "partitionfuns.h"
#include "distances.h"
#include "misc.h"
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


void update_dist_mat(mat& dist_mat, const mat& rankings,
                     const mat& rho_old, const std::string& metric,
                     const vec& obs_freq){
  int n_clusters = dist_mat.n_cols;
  for(int i = 0; i < n_clusters; ++i)
    dist_mat.col(i) = rank_dist_vec(rankings, rho_old.col(i), metric, obs_freq);
}

uvec update_cluster_labels(
    const mat& dist_mat,
    const vec& cluster_probs,
    const vec& alpha_old,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<vec> logz_estimate = R_NilValue,
    const bool& save_ind_clus = false
){
  int n_assessors = dist_mat.n_rows;
  int n_clusters = dist_mat.n_cols;
  uvec new_cluster_assignment(n_assessors);


  mat assignment_prob(n_assessors, n_clusters);
  for(int i = 0; i < n_clusters; ++i){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      alpha_old(i) / n_items * dist_mat.col(i) -
      get_partition_function(n_items, alpha_old(i), cardinalities, logz_estimate, metric);
  }

  for(int i = 0; i < n_assessors; ++i){
    // Exponentiate to get unnormalized prob relative to max
    rowvec probs = exp(assignment_prob.row(i) -
      max(assignment_prob.row(i)));

    // Normalize with 1-norm
    probs = arma::normalise(probs, 1);

    assignment_prob.row(i) = probs;
    new_cluster_assignment(i) = sample_int(assignment_prob.row(i));
  }

  if(save_ind_clus){
    assignment_prob.save(std::string("cluster_probs") + std::to_string(t + 1) + std::string(".csv"), arma::csv_ascii);
  }
  return(new_cluster_assignment);
}

vec update_cluster_probs(
    const uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
){

  vec cluster_probs(n_clusters);

  for(int i = 0; i < n_clusters; ++i){
    // Find the parameter for this cluster and provide it to the gamma distribution
    cluster_probs(i) = R::rgamma(sum(current_cluster_assignment == i) + psi, 1.0);
  }
  // Finally, normalize cluster_probs with 1-norm.
  // result now comes from Dirichlet(tau_k(0), ..., tau_k(n_clusters))
  return arma::normalise(cluster_probs, 1);

}


vec update_wcd(const uvec& current_cluster_assignment,
                     const mat& dist_mat){
  int n_clusters = dist_mat.n_cols;
  vec wcd(n_clusters);

  uvec inds = arma::regspace<uvec>(0, n_clusters - 1);
  for(int i = 0; i < n_clusters; ++i){
    mat dist_vec = dist_mat.submat(arma::find(current_cluster_assignment == i), inds.subvec(i, i));
    wcd(i) = sum(conv_to<vec>::from(dist_vec));
  }

  return(wcd);
}
