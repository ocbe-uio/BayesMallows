#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec setdiff(const arma::vec& v1, const arma::vec& v2) {
  Rcpp::NumericVector v = Rcpp::setdiff(Rcpp::as<Rcpp::NumericVector>(v1), Rcpp::as<Rcpp::NumericVector>(v2));
  return v;
}

/*** R
v1 <- 1:10
v2 <- 1:3
setdiff(v1, v2)
*/
