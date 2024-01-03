#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
cube abind(cube x, cube y) {
  return join_slices(x, y);
}

/*** R
a <- array(2, dim = c(2, 2, 4))
b <- array(3, dim = c(2, 2, 2))
abind(a, b)
*/
