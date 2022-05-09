// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

arma::vec sample(arma::vec x, int size, bool replace = false){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}

arma::uvec sample(arma::uvec x, int size, bool replace = false){
  return Rcpp::RcppArmadillo::sample(x, size, replace);
}
