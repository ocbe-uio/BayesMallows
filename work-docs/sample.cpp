#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
uvec foo() {
  vec a = Rcpp::runif(10, 0, 1);
  Rcpp::Rcout << "a = " << a << std::endl;
  uvec b = sort_index(Rcpp::as<vec>(Rcpp::wrap(a)));
  return b;
}

/*** R
set.seed(1)
foo()
*/
