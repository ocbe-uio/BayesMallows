#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_alpha(arma::vec&, arma::vec& ,
                  double&,
                  const arma::mat&, int& ,
                  const arma::vec& ,
                  const double&, const std::string&,
                  const double&, const int&, const int&,
                  Rcpp::Nullable<arma::vec>,
                  Rcpp::Nullable<arma::vec>);

void update_rho(arma::mat& rho, arma::vec& rho_acceptance, arma::vec& rho_old,
                const double& alpha_old, const int& L, const arma::mat& R,
                const std::string& metric, const int& n, const int& t);

#endif
