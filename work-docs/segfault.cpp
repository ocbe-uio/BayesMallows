#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
double takeLog3(double val) {
  if (val <= 0.0) {           // log() not defined here
    stop("Inadmissible value");
  }
  return log(val);
}

/*** R
takeLog3(-1)
*/
