#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"

arma::mat initialize_rho(int n_items, int n_cols,
                         Rcpp::Nullable<arma::mat> rho_init = R_NilValue);
double update_alpha(arma::vec& alpha_acceptance,
                  const double& alpha_old,
                  const arma::mat& rankings,
                  const arma::vec& obs_freq,
                  const int& cluster_index,
                  const arma::vec& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                  double alpha_max = 1e6);

void update_rho(arma::cube& rho, arma::vec& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& rho_thinning,
                const double& alpha_old, const int& leap_size, const arma::mat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices, const arma::vec& obs_freq);



#endif
