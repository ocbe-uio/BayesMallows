#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec normalize_weights(const vec& log_inc_wgt){
  double maxw = max(log_inc_wgt);
  double log_sum_exp = maxw + log(sum(exp(log_inc_wgt - maxw)));
  return exp(log_inc_wgt - log_sum_exp);
}

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
  return(prob / tgamma(set_length + 1));
}

bool is_pseudo(const std::string aug_method, const std::string metric) {
  // Checks for valid combinations of the inputs, stops if invalid
  if (aug_method == "uniform") {
    return(false);
  } else if (aug_method == "pseudo") {
    if ((metric == "footrule") || (metric == "spearman")) {
      return(true);
    } else {
      Rcpp::stop("Pseudolikelihood only supports footrule and spearman metrics");
    }
  } else {
    Rcpp::stop("Invalid aug_method. Please choose random or pseudolikelihood");
  }
}
