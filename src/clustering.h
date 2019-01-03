#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "RcppArmadillo.h"
#include "partitionfuns.h"
#include "misc.h"

void update_cluster_labels(
    arma::uvec& current_cluster_assignment,
    const arma::mat& dist_mat,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue
);

arma::vec update_cluster_probs(
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
);

arma::vec update_wcd(const arma::uvec& current_cluster_assignment,
                     const arma::mat& dist_mat);

#endif
