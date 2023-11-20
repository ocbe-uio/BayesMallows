#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"

arma::mat initialize_rho(int n_items, int n_cols,
                         Rcpp::Nullable<arma::mat> rho_init = R_NilValue);
double update_alpha(
                  const double& alpha_old,
                  const arma::mat& rankings,
                  const arma::vec& observation_frequency,
                  const arma::vec& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::List& logz_list);

arma::vec make_new_rho(arma::vec current_rho, const arma::mat& rankings, double alpha_old, int leap_size, std::string metric,
                 arma::vec observation_frequency);

#endif
