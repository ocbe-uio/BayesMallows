#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int test(Rcpp::List l) {
  Rcpp::IntegerVector a = Rcpp::as<Rcpp::IntegerVector>(l["a"]);


  if(IntegerVector::is_na(a[0])) {
    Rcpp::Rcout << "hei" << std::endl;
  }
  return 0;
}

/*** R
test(list(a = c(NA, 1)))
*/
