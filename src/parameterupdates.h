#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"

arma::umat initialize_rho(uint n_items, uint n_cols,
                         Rcpp::Nullable<arma::umat> rho_init = R_NilValue);
double update_alpha(arma::vec& alpha_acceptance,
                  const double& alpha_old,
                  const arma::umat& rankings,
                  const arma::uvec& obs_freq,
                  const int& cluster_index,
                  const arma::uvec& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                  double alpha_max = 1e6);

void update_rho(arma::ucube& rho, arma::vec& rho_acceptance, arma::umat& rho_old,
                int& rho_index, const int& cluster_index, const int& rho_thinning,
                const double& alpha_old, const int& leap_size, const arma::umat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices, const arma::uvec& obs_freq);



#endif
