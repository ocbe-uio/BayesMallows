#ifndef SMC_H
#define SMC_H

#include "RcppArmadillo.h"

arma::vec normalize_weights(const arma::vec& log_inc_wgt);
arma::vec initialize_alpha(const int& n_particles, const Rcpp::Nullable<arma::vec>& alpha_init = R_NilValue);
double get_exponent_sum(double, arma::vec, int, arma::mat, std::string);

#endif
