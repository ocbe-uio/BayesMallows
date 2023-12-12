#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerVector test(int size) {
  return Rcpp::sample(10, 1);
}

/*** R
set.seed(1)
test(100)
#apply(sapply(1:100, function(i) test(100)), 1, mean)
*/
