#include <RcppArmadillo.h>
#include <Rmath.h>
using namespace arma;

// [[Rcpp::export]]
uvec rmn(vec probs) {
  int k = probs.size();
  ivec ans(k);
  R::rmultinom(1, probs.begin(), k, ans.begin());
  uvec ret = find(ans == 1);
  return(ret);
}


/*** R
rmn(c(.005, 0, .005, .99))
*/
