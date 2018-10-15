#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>

int factorial(int);
int binomial_coefficient(int, int);
arma::uvec std_setdiff(arma::uvec&, arma::uvec&);
int sample_int(const arma::rowvec& probs);

#endif
