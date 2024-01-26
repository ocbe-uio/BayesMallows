#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ouch() {
  NumericVector x(10);
  x[1000000] = 1;
  return x;
}

/*** R
ouch()
*/
