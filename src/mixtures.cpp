#include "RcppArmadillo.h"
#include <cmath>
#include "partitionfuns.h"
#include <string>

// [[Rcpp::depends(RcppArmadillo)]]


void update_dist_mat(arma::mat& dist_mat, const arma::mat& rankings, const arma::mat& rho_old, const std::string& metric){
  int n_clusters = dist_mat.n_cols;
  for(int i = 0; i < n_clusters; ++i)
    dist_mat.col(i) = rank_dist_vec(rankings, rho_old.col(i), metric);
}

arma::uvec update_cluster_labels(
    const arma::mat& dist_mat,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
    const bool& save_ind_clus = false
){
  int n_assessors = dist_mat.n_rows;
  int n_clusters = dist_mat.n_cols;
  arma::uvec new_cluster_assignment(n_assessors);


  arma::mat assignment_prob(n_assessors, n_clusters);
  for(int i = 0; i < n_clusters; ++i){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      alpha_old(i) / n_items * dist_mat.col(i) -
      get_partition_function(n_items, alpha_old(i), cardinalities, logz_estimate, metric);
  }

  for(int i = 0; i < n_assessors; ++i){
    // Exponentiate to get unnormalized prob relative to max
    arma::rowvec probs = arma::exp(assignment_prob.row(i) -
      arma::max(assignment_prob.row(i)));

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

arma::vec update_cluster_probs(
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
){

  arma::vec cluster_probs(n_clusters);

  for(int i = 0; i < n_clusters; ++i){
    // Find the parameter for this cluster and provide it to the gamma distribution
    cluster_probs(i) = R::rgamma(arma::sum(current_cluster_assignment == i) + psi, 1.0);
  }
  // Finally, normalize cluster_probs with 1-norm.
  // result now comes from Dirichlet(tau_k(0), ..., tau_k(n_clusters))
  return arma::normalise(cluster_probs, 1);

}


arma::vec update_wcd(const arma::uvec& current_cluster_assignment,
                     const arma::mat& dist_mat){
  int n_clusters = dist_mat.n_cols;
  arma::vec wcd(n_clusters);

  arma::uvec inds = arma::regspace<arma::uvec>(0, n_clusters - 1);
  for(int i = 0; i < n_clusters; ++i){
    arma::mat dist_vec = dist_mat.submat(arma::find(current_cluster_assignment == i), inds.subvec(i, i));
    wcd(i) = arma::sum(arma::conv_to<arma::vec>::from(dist_vec));
  }

  return(wcd);
}
