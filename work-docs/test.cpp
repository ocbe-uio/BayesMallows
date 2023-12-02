#include <Rcpp.h>
using namespace Rcpp;

struct Data {
  Data(List data) :
  constraints { as<List>(data["constraints"]) } {}
  ~Data() = default;

  const List constraints;
};

// [[Rcpp::export]]
int foo(List dat) {
  Data dd {dat};
  int result = dd.constraints.length();
  return result;
}

/*** R
data <- list(a = rnorm(100), constraints = NULL, b = letters)
foo(data)
*/
