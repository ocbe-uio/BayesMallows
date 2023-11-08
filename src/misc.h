#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h> // need because of functions that return arma and Rcpp objects

double rtruncbeta(int shape1, int shape2, double trunc = 1);
double divide_by_fact(double, int);
arma::vec normalize_weights(const arma::vec& log_inc_wgt);

#endif
