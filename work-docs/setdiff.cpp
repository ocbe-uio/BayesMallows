#include <Rcpp.h>

// [[Rcpp::export]]
void hello(Rcpp::CharacterVector v1, Rcpp::CharacterVector v2) {
  Rcpp::Rcout << Rcpp::setdiff(v1, v2) << std::endl;
}

/*** R
hello(c("hello", "hei", "a", "b"), c("b", "c"))
*/
