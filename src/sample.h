#ifndef SAMPLE_H
#define SAMPLE_H

#include <RcppArmadillo.h>
arma::vec sample(arma::vec x, int size, bool replace = false);
arma::uvec sample(arma::uvec x, int size, bool replace = false);

#endif
