#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector createEmptyIntegerVector() {
  IntegerVector x = IntegerVector();
  for(auto i : x) {Rcpp::Rcout << i << std::endl;}
    return IntegerVector(); // Creates an empty IntegerVector
}

/*** R
# Example usage
empty_vector <- createEmptyIntegerVector()
print(empty_vector)
*/
