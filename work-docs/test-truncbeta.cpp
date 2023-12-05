#include <Rcpp.h>
#include <Rmath.h>
using namespace Rcpp;

//[[Rcpp::export]]
double rtruncbeta(int shape1, int shape2, double trunc = 1) {
  for(size_t i{}; i < 1e8; i++){
    double x{R::rbeta(shape1, shape2)};
    if(x < trunc) return x;
  }
  Rcpp::stop("Unable to sample from truncated beta distribution.\n");
}


/*** R
sample <- replicate(10000, rtruncbeta(1, 3, .5))
hist(sample)
*/
