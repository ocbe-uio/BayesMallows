#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>

void augment_pairwise(
    arma::mat& rankings,
    const arma::umat& cluster_indicator,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const arma::mat& pairwise_preferences,
    const arma::mat& constrained_elements,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::mat& aug_acceptance,
    int& aug_diag_index,
    const int& aug_diag_thinning
);

#endif
