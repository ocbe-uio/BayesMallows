#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int test(List x) {
  for(Rcpp::List y : x) {
    std::vector<double> z(y["number"]);
    for(auto w : z) Rcpp::Rcout << w << std::endl;
  }
  return 0;
}

/*** R
test(list(list(number = 1, cumber = "some"), list(cumber = "cu", number = rnorm(10))))
test(list())
*/
