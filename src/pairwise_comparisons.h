#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

void update_shape_bernoulli(
    double& shape_1,
    double& shape_2,
    const double& kappa_1,
    const double& kappa_2,
    const arma::mat& rankings,
    const Rcpp::List& constraints
);

void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const double& theta,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    arma::vec& aug_acceptance,
    const bool& clustering,
    const std::string& error_model,
    const int& Lswap
);


#endif
