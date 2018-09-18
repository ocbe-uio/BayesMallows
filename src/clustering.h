#ifndef CLUSTERING_H
#define CLUSTERING_H

#include "RcppArmadillo.h"
#include <math.h>
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_cluster_labels(
    arma::umat& cluster_assignment,
    arma::uvec& current_cluster_assignment,
    const arma::mat& dist_mat,
    const arma::mat& rho_old,
    const arma::mat& rankings,
    const arma::mat& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const int& n_assessors,
    const int& n_clusters,
    const int& cluster_assignment_thinning,
    int& cluster_assignment_index,
    const int& t,
    const std::string& metric,
    const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
    const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue
);

void update_cluster_probs(
        arma::mat& cluster_probs,
        const arma::uvec& current_cluster_assignment,
        const int& n_clusters,
        const int& psi,
        const int& t
);

void update_wcd(arma::mat& within_cluster_distance, const arma::uvec& current_cluster_assignment,
                const arma::mat& dist_mat, const int& n_clusters, const int& t);

#endif
