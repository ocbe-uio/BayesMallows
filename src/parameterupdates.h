#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_alpha(arma::mat& alpha,
                  arma::vec& alpha_acceptance,
                  arma::vec& alpha_old,
                  const arma::mat& rankings,
                  int& alpha_index,
                  int& cluster_index,
                  const arma::mat& rho_old,
                  const double& sd_alpha,
                  const std::string& metric,
                  const double& lambda,
                  const int& n_items,
                  const Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                  const Rcpp::Nullable<arma::vec> is_fit = R_NilValue);

void update_rho(arma::cube& rho, arma::vec& rho_acceptance, arma::mat& rho_old,
                int& rho_index, const int& cluster_index, const int& thinning,
                const double& alpha_old, const int& leap_size, const arma::mat& rankings,
                const std::string& metric, const int& n_items, const int& t,
                const arma::uvec& element_indices, bool& rho_accepted);


#endif
