#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::cube abind(arma::cube x, arma::cube y) {
  return arma::join_slices(x, y);
}
