#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec test() {
  arma::vec probs(20);
  probs.fill(.05);
  arma::uvec inds = arma::regspace<arma::uvec>(0, 19);
  return Rcpp::RcppArmadillo::sample(inds, inds.size(), true, probs);
}


/*** R
set.seed(3)
test()
*/
