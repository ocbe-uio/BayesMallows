#ifndef MISC_H
#define MISC_H

#include <RcppArmadillo.h> // need because of functions that return arma and Rcpp objects
using namespace arma;

long int factorial(int);
int binomial_coefficient(int, int);
uvec std_setdiff(uvec&, uvec&);
int sample_int(const rowvec& probs);
double rtruncbeta(int shape1, int shape2, double trunc = 1);
uvec arma_setdiff(uvec x, uvec y);
vec arma_setdiff_vec(vec, vec, const bool& = false);
Rcpp::NumericVector Rcpp_setdiff_arma(ivec, vec);
uvec maybe_offset_indices(vec&, uvec, const bool& = true);
sword sample_one_with_prob(vec, vec);
uvec new_pseudo_proposal(uvec);
double divide_by_fact(double, int);
uvec permute_with_weights(vec, int);
vec arma_vec_seq(int);
#endif
