#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distfuns.h"

void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::vec& aug_acceptance,
    const bool& clustering,
    bool& augmentation_accepted
);


#endif
