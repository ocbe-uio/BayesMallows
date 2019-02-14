#ifndef MISSING_H
#define MISSING_H

#include <RcppArmadillo.h>
#include "distances.h"
#include "misc.h"

void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator,
                              const arma::uvec& assessor_missing);

void update_missing_ranks(arma::mat& rankings, const arma::uvec& current_cluster_assignment,
                          arma::vec& aug_acceptance,
                          const arma::umat& missing_indicator,
                          const arma::uvec& assessor_missing,
                          const arma::vec& alpha, const arma::mat& rho,
                          const std::string& metric);

#endif
