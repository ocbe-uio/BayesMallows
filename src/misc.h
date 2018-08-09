#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h>
#include <algorithm>

int factorial(int);
int binomial_coefficient(int, int);
arma::uvec std_setdiff(arma::uvec&, arma::uvec&);

#endif
