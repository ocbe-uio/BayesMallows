#ifndef LEAPANDSHIFT_H
#define LEAPANDSHIFT_H

#include "RcppArmadillo.h"

Rcpp::List leap_and_shift(const arma::vec&, int);
void shift_step(arma::vec& proposal, const arma::vec& rho,
                const int& u, double& delta_r, arma::vec& indices);

#endif
