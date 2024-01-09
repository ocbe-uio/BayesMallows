#include <RcppArmadillo.h>
#include <execution>
using namespace arma;

// [[Rcpp::export]]
vec setdiff(const vec& x, const vec& y) {
  vec xs = sort(x);
  vec ys = sort(y);

  std::vector<double> diff;
  std::set_difference(
    std::execution::par,
    xs.begin(), xs.end(), ys.begin(), ys.end(),
    std::inserter(diff, diff.begin()));

  return conv_to<vec>::from(diff);
}

// [[Rcpp::export]]
vec setdiff2(const vec& x, const vec& y) {
  vec xs = sort(x);
  vec ys = sort(y);

  std::vector<double> diff;
  std::set_difference(
    std::execution::par_unseq,
    xs.begin(), xs.end(), ys.begin(), ys.end(),
    std::inserter(diff, diff.begin()));

  return conv_to<vec>::from(diff);
}

/*** R
microbenchmark::microbenchmark(
  setdiff(sample(100, 100), sample(100, 10)),
  setdiff2(sample(100, 100), sample(100, 10))
)
*/
