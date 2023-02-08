#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h> // need because of functions that return arma and Rcpp objects

double rtruncbeta(int shape1, int shape2, double trunc = 1);
arma::uvec maybe_offset_indices(arma::vec&, arma::uvec, const bool& = true);
double divide_by_fact(double, int);
bool is_pseudo(const std::string, const std::string);
#endif
