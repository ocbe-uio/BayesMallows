#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int test() {
  return arma::randi<int>(arma::distr_param(1, 20));
}

/*** R
set.seed(1)
mean(replicate(1e6, test()))
*/
