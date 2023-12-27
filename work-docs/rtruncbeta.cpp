#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rtruncbeta1(int shape1, int shape2, double trunc = .5) {
  for(size_t i{}; i < 1e3; i++){
    double x{R::rbeta(shape1, shape2)};
    if(x < trunc) return x;
  }
  Rcpp::stop("Unable to sample from truncated beta distribution.\n");
}

// [[Rcpp::export]]
double rtruncbeta2(double shape1, double shape2, double trunc = .5) {
  double Fa = R::pbeta(0, shape1, shape2, true, false);
  double Fb = R::pbeta(trunc, shape1, shape2, true, false);
  double u = R::runif(Fa, Fb);
  return R::qbeta(u, shape1, shape2, true, false);
}


/*** R

v1 <- replicate(1000, rtruncbeta1(1, 2))
v2 <- replicate(1000, rtruncbeta2(1, 2))
*/
