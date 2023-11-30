#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
std::vector<double> setdiff(const arma::vec& v1, const arma::vec& v2) {
  std::vector<double> diff;
  arma::vec s1 = sort(v1);
  arma::vec s2 = sort(v2);
  std::set_difference(s1.begin(), s1.end(),
                      s2.begin(), s2.end(),
                      std::inserter(diff, diff.begin()));

  return diff;
}

/*** R
v1 <- sample(1:10)
v2 <- sample(1:3)
setdiff(v1, v2)
*/
