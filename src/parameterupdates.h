#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "misc.h"
#include "leapandshift.h"
#include "distfuns.h"
#include "partitionfuns.h"

void update_alpha(arma::vec&,
                  arma::vec&,
                  double&,
                  const arma::mat&,
                  int&,
                  const arma::vec&,
                  const double&,
                  const std::string&,
                  const double&,
                  const int&,
                  const int&,
                  Rcpp::Nullable<arma::vec>,
                  Rcpp::Nullable<arma::vec>);

void update_rho(arma::mat&,
                arma::vec&,
                arma::vec&,
                int&,
                const int&,
                const double&,
                const int&,
                const arma::mat&,
                const std::string&,
                const int&,
                const int&);

#endif
