#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Truncated beta distribution
double rtruncbeta(int shape1, int shape2, double trunc = 1) {
  int i = 0;
  double x;
  while(i < 1000){
    x = chi2rnd(2 * shape1);
    x = x / (x + chi2rnd(2 * shape2));

    if(x < trunc) break;
    ++i;
  }
  return x;
}

uvec maybe_offset_indices(
  vec& x,
  uvec idx_x,
  const bool& quiet = true
) {
  // Adjust the indices of x (i.e., idx_x) depending on whether it seems to be
  // using R or C++ indices.
  const uvec& io_idx_cpp   = find_nonfinite(x);
  const uvec& io_idx_input = sort(idx_x);
  std::string message = "C++ indices detected. Unchanged.";
  if (any(io_idx_input - io_idx_cpp)) {
    idx_x -= 1;
    message = "R indices detected. Shifted.";
  }
  if (!quiet) {
    Rcpp::Rcout << message << std::endl;
  }
  return(idx_x);
}


double divide_by_fact(double prob, int set_length) {
  // Using the fact that Gamma(x + 1) = x!
  prob /= tgamma(set_length + 1);
  return(prob);
}
