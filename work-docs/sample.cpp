#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
uvec foo() {
  ivec a1 = regspace<ivec>(0, 3);
  ivec a = Rcpp::sample(a1.size(), 2);
  Rcpp::Rcout << "a = " << a << std::endl;
  uvec b = sort_index(Rcpp::as<vec>(Rcpp::wrap(a)));
  return b;
}

/*** R
set.seed(1)
foo()
*/
