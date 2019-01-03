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

void update_cluster_probs(
        arma::mat& cluster_probs,
        const arma::uvec& current_cluster_assignment,
        const int& psi,
        const int& t
);

void update_wcd(arma::mat& within_cluster_distance, const arma::uvec& current_cluster_assignment,
                const arma::mat& dist_mat, const int& n_clusters, const int& t);

#endif
