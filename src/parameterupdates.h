#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_alpha(arma::mat& alpha,
                  arma::vec& alpha_acceptance,
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
                  Rcpp::Nullable<arma::vec> is_fit = R_NilValue);

void update_rho(arma::cube& rho, arma::vec& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& thinning,
                const double& alpha_old, const int& L, const arma::mat& R,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices);

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
);

void update_cluster_probs(
        arma::mat& cluster_probs,
        const arma::umat& cluster_indicator,
        const int& n_clusters,
        const int& psi,
        const int& t
);

#endif
