#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
arma::ivec match(arma::ivec v1, arma::ivec v2) {
  ivec indices(v1.n_elem);

  // Loop through each element of v1
  for (int i = 0; i < v1.n_elem; ++i) {
    // Find the index of v1[i] in v2
    indices(i) = index(v2, v1(i));
  }

  return indices;
}

/*** R
test(list(b = 3))
*/
