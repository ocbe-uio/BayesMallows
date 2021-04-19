#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

double get_mallows_loglik(double, arma::vec, int, arma::mat, std::string);
Rcpp::List leap_and_shift_probs(arma::vec rho, int leap_size, int n_items);
arma::vec metropolis_hastings_rho(double alpha, int n_items, arma::mat rankings, std::string metric, arma::vec rho, int leap_size);
double metropolis_hastings_alpha(double alpha, int n_items, arma::mat rankings, std::string metric, arma::vec rho, const Rcpp::Nullable<arma::vec> logz_estimate, double alpha_prop_sd, double lambda, double alpha_max);
#endif
