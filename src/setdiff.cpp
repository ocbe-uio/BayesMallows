#include "RcppArmadillo.h"
using namespace arma;

vec setdiff(const vec& x, const vec& y) {
  vec xs = sort(x);
  vec ys = sort(y);

  std::vector<double> diff;
  std::set_difference(
    xs.begin(), xs.end(), ys.begin(), ys.end(),
    std::inserter(diff, diff.begin()));

  return conv_to<vec>::from(diff);
}
