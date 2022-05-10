#ifndef SAMPLE_H
#define SAMPLE_H

#include <RcppArmadillo.h>
arma::vec sample(const arma::vec& x, const int& size, const bool& replace, const arma::vec& probs);
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace, const arma::vec& probs);
#endif
