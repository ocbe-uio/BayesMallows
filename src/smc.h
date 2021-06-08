#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

double get_mallows_loglik(double, arma::vec, int, arma::mat, std::string);
Rcpp::List leap_and_shift_probs(arma::vec, int, int);
arma::vec metropolis_hastings_rho(double, int, arma::mat, std::string, arma::vec, int);
double metropolis_hastings_alpha(double, int, arma::mat, std::string, arma::vec, const Rcpp::Nullable<arma::vec>, double, double, double);

#endif
