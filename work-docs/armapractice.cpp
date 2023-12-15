#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins("cpp11")]]

using namespace arma;
// [[Rcpp::export]]
arma::vec colCumSum(arma::mat& X) {
  X.each_col( [](vec& a){ double b = sum(a); } );
  return X.col(0);
}

/*** R
M <- matrix(1:16, 4, 4)
colCumSum(M)
*/
