#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>

void augment_pairwise(
    arma::mat& R,
    const double& alpha,
    const arma::vec& rho,
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
