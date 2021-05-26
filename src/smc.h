#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

double get_mallows_loglik(double, arma::vec, int, arma::mat, std::string);
Rcpp::List leap_and_shift_probs(arma::vec rho, int leap_size, int n_items);

#endif
