#ifndef SAMPLE_H
#define SAMPLE_H

#include <RcppArmadillo.h>
arma::vec sample(arma::vec x, int size, bool replace = false);
arma::vec sample(arma::vec x, int size, bool replace, arma::vec probs);
arma::uvec sample(arma::uvec x, int size, bool replace = false);
arma::uvec sample(arma::uvec x, int size, bool replace, arma::vec probs);


#endif
