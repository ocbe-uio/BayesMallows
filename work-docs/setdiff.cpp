#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int test(Rcpp::List a)  {
  auto tmp = Rcpp::IntegerVector(a["b"]);
  int ret;
  if(tmp[0] == R_NilValue) {
    ret = 0;
  } else {
    ret = tmp[0];
  }
  return ret;
}

/*** R
test(list(b = 3))
*/
