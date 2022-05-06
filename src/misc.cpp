#include <RcppArmadillo.h>
#include <cmath>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to sample an integer given a set of probabilities
// [[Rcpp::export]]
int sample_int(const arma::rowvec& probs){

  if(probs.has_nan() || probs.has_inf()){
    Rcpp::Rcout << "probs = " << probs << std::endl;
    Rcpp::stop("Cannot sample_int.");
  }

  // Draw a uniform random number
  double u = randu();

  uvec matches = find(cumsum(probs) > u, 1, "first");

  return as_scalar(matches);
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

// This is practically a signed variation of arma_setdiff. It also uses Rcpp as
// a crutch, so future changes that work around this and only use arma objects
// are welcome.
vec arma_setdiff_vec(vec x, vec y, const bool& sort_unique = false) {
  Rcpp::NumericVector x_Rcpp, y_Rcpp;
  vec x_y_diff;
  x_Rcpp = x;
  y_Rcpp = y;
  if (sort_unique) {
    x_y_diff = Rcpp::sort_unique(Rcpp::setdiff(x_Rcpp, y_Rcpp));
  } else {
    x_y_diff = Rcpp::setdiff(x_Rcpp, y_Rcpp);
  }
  return x_y_diff;
}

// This is practically an Rcpp variation of arma_setdiff_arma. Future changes
// that eliminate the use of this function in favor of arma_setdiff_vec() are
// welcome.
Rcpp::NumericVector Rcpp_setdiff_arma(ivec x, vec y) {
  Rcpp::NumericVector x_Rcpp, y_Rcpp;
  x_Rcpp = x;
  y_Rcpp = y;
  Rcpp::NumericVector x_y_diff = Rcpp::setdiff(x_Rcpp, y_Rcpp);
  return x_y_diff;
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

sword sample_one_with_prob(vec set, vec probs) {
  // Used in SMC to fill in the new augmented ranking going forward.
  Rcpp::NumericVector set_Rcpp, probs_Rcpp;
  set_Rcpp = set;
  probs_Rcpp = probs;
  sword chosen_one = Rcpp::as<int>(Rcpp::sample(set_Rcpp, 1, false, probs_Rcpp));
  return(chosen_one);
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

