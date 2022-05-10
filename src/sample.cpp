// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>
arma::vec sample(const arma::vec& x, const int& size, const bool& replace, const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}
arma::uvec sample(const arma::uvec& x, const int& size, const bool& replace, const arma::vec& probs){
  return Rcpp::RcppArmadillo::sample(x, size, replace, probs);
}
