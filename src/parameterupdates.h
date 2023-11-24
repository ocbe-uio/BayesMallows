#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"

arma::mat initialize_rho(int n_items, int n_cols,
                         Rcpp::Nullable<arma::mat> rho_init = R_NilValue);

arma::vec make_new_rho(arma::vec current_rho, const arma::mat& rankings, double alpha_old, int leap_size, std::string metric,
                 arma::vec observation_frequency);

#endif
