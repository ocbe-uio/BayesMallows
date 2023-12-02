#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int list_of_list_length(List x) {
  List y = as<List>(x["a"]);
  return y.length();
}

/*** R
ll <- list()
list_of_list_length(ll)

*/
