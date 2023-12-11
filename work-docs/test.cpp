#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec test() {
  return arma::randu(1);
}


/*** R
set.seed(1)
test()
*/
