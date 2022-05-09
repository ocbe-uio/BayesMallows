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

uvec new_pseudo_proposal(uvec items) {
  // Used in SMC to create new agumented ranking by using pseudo proposal. This
  // function randomly permutes the unranked items to give the order in which
  // they will be allocated.
  Rcpp::IntegerVector items_Rcpp;
  items_Rcpp = items;
  uvec order;
  order = Rcpp::as<uvec>(Rcpp::sample(items_Rcpp, items_Rcpp.length())) + 1;
  return(order);
}

double divide_by_fact(double prob, int set_length) {
  // Using the fact that Gamma(x + 1) = x!
  prob /= tgamma(set_length + 1);
  return(prob);
}

uvec permute_with_weights(vec weights, int N) {
  // Using weights_Rcpp so that Rcpp::sample compiles. More details on
  // https://github.com/ocbe-uio/BayesMallows/issues/90#issuecomment-866614296
  Rcpp::NumericVector weights_Rcpp;
  weights_Rcpp = weights;
  uvec index;
  index = Rcpp::as<uvec>(Rcpp::sample(N, N, true, weights_Rcpp)) - 1;
  return(index);
}

