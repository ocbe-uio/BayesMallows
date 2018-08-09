#ifndef MISSING_H
#define MISSING_H

#include <RcppArmadillo.h>

void define_missingness(arma::mat& missing_indicator, arma::vec& assessor_missing,
                        const arma::mat& R,
                        const int& n_items, const int& n_assessors);

arma::vec propose_augmentation(const arma::vec& ranks, const arma::vec& indicator,
                               const int& n_items);

void initialize_missing_ranks(arma::mat& R, const arma::mat& missing_indicator,
                              const arma::vec& assessor_missing,
                              const int& n_items, const int& n_assessors);

void update_missing_ranks(arma::mat& R, arma::mat& aug_acceptance,
                          const arma::mat& missing_indicator,
                          const arma::vec& assessor_missing,
                          const int& n_items, const int& n_assessors,
                          const double& alpha, const arma::vec& rho,
                          const std::string& metric, const int& t);

#endif
