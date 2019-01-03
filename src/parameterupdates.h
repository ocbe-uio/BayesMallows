#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "leapandshift.h"
#include "distances.h"
#include "partitionfuns.h"

double rtruncbeta(int shape1, int shape2, double trunc = 1);

void update_alpha(arma::mat& alpha,
                  arma::vec& alpha_acceptance,
                  arma::vec& alpha_old,
                  const arma::mat& rankings,
                  int& alpha_index,
                  int& cluster_index,
                  const arma::mat& rho_old,
                  const double& alpha_prop_sd,
                  const std::string& metric,
                  const double& lambda,
                  const int& n_items,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue);

void update_rho(arma::cube& rho, arma::vec& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& rho_thinning,
                const double& alpha_old, const int& leap_size, const arma::mat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices, bool& rho_accepted);


#endif
