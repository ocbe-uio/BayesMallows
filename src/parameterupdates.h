#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_alpha(arma::vec& alpha, arma::vec& alpha_acceptance,
                  double& alpha_old,
                  const arma::mat& R, int& alpha_index,
                  const arma::vec& rho_old,
                  const double& sd_alpha, const std::string& metric,
                  const double& lambda, const int& n, const int& N,
                  const arma::vec& cardinalities);

void update_rho(arma::mat& rho, arma::vec& rho_acceptance, arma::vec& rho_old,
                const double& alpha_old, const int& L, const arma::mat& R,
                const std::string& metric, const int& n, const int& t);

#endif
