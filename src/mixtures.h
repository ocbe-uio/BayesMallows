#ifndef MIXTURES_H
#define MIXTURES_H

#include "RcppArmadillo.h"

void update_dist_mat(arma::mat& dist_mat, const arma::mat& rankings,
                     const arma::mat& rho_old, const std::string& metric,
                     const arma::vec& observation_frequency);

arma::uvec update_cluster_labels(
    const arma::mat& dist_mat,
    const arma::vec& cluster_probs,
    const arma::vec& alpha_old,
    const int& n_items,
    const int& t,
    const std::string& metric,
    const Rcpp::List& logz_list,
    const bool& save_ind_clus = false
);

arma::vec update_cluster_probs(
    const arma::uvec& current_cluster_assignment,
    const int& n_clusters,
    const int& psi
);

arma::vec update_wcd(const arma::uvec& current_cluster_assignment,
                     const arma::mat& dist_mat);

#endif
