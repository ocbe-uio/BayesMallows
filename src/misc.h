#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h> // need because of functions that return arma and Rcpp objects

arma::uvec std_setdiff(arma::uvec&, arma::uvec&);
int sample_int(const arma::rowvec& probs);
double rtruncbeta(int shape1, int shape2, double trunc = 1);
arma::uvec arma_setdiff(arma::uvec x, arma::uvec y);
arma::vec arma_setdiff_vec(arma::vec, arma::vec, const bool& = false);
Rcpp::NumericVector Rcpp_setdiff_arma(arma::ivec, arma::vec);
arma::uvec maybe_offset_indices(arma::vec&, arma::uvec, const bool& = true);
arma::sword sample_one_with_prob(arma::vec, arma::vec);
arma::uvec new_pseudo_proposal(arma::uvec);
double divide_by_fact(double, int);
arma::uvec permute_with_weights(arma::vec, int);
arma::vec arma_vec_seq(int);
#endif
