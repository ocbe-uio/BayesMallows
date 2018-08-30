#ifndef PAIRWISE_H
#define PAIRWISE_H

#include <RcppArmadillo.h>

void augment_pairwise(
    arma::mat& rankings,
    const arma::umat& cluster_assignment,
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

void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const arma::vec& current_ranking);

void propose_pairwise_augmentation(arma::vec& proposal,
                                   const arma::mat& rankings,
                                   const Rcpp::List& constraints,
                                   const int& n_items,
                                   const int& i);


#endif
