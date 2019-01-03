#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

void update_shape_bernoulli(
    arma::vec& shape_1,
    arma::vec& shape_2,
    const double& kappa_1,
    const double& kappa_2,
    const int& n_assessors,
    const int& n_items,
    const arma::mat& rankings,
    const Rcpp::List& constraints,
    const int& t
);

void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const double& theta,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    const int& n_assessors,
    const int& n_items,
    const int& t,
    arma::vec& aug_acceptance,
    const bool& clustering,
    bool& augmentation_accepted,
    std::string error_model
);


#endif
