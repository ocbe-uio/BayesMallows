#include <Rcpp.h>
#include <Rmath.h>
using namespace std;

// [[Rcpp::export]]
double s() {
  return R::unif_rand();
}


/*** R
s()
*/
