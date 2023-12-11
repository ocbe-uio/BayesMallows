#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec test() {
  arma::vec probs(5);
  probs.fill(.2);
  arma::uvec inds = arma::regspace<arma::uvec>(0, 4);
  return Rcpp::RcppArmadillo::sample(inds, inds.size(), true, probs);
}


/*** R
set.seed(1)
test()
*/
