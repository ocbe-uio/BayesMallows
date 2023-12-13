#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector test(int size) {
  return sample(size, size) - 1;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
set.seed(1)
test(10)
*/
