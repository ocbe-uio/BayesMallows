#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function to compute the factorial
// taken from http://www.cplusplus.com/forum/unices/33379/
// [[Rcpp::export]]
long int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

// [[Rcpp::export]]
int binomial_coefficient(int n, int k){

  // Special case:
  if( k > n ) return 0;

  int res = 1;

  // Since C(n, k) = C(n, n-k)
  if ( k > n - k )
    k = n - k;

  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
  for (int i = 0; i < k; ++i)
  {
    res *= (n - i);
    res /= (i + 1);
  }

  return res;
}


// Function to sample an integer given a set of probabilities
// [[Rcpp::export]]
int sample_int(const arma::rowvec& probs){

  if(probs.has_nan() || probs.has_inf()){
    Rcpp::Rcout << "probs = " << probs << std::endl;
    Rcpp::stop("Cannot sample_int.");
  }

  // Draw a uniform random number
  double u = randu();

  uvec matches = find(arma::cumsum(probs) > u, 1, "first");

  return arma::as_scalar(matches);
}

// Truncated beta distribution
double rtruncbeta(int shape1, int shape2, double trunc = 1) {
  int i = 0;
  double x;
  while(i < 1000){
    x = arma::chi2rnd(2 * shape1);
    x = x / (x + arma::chi2rnd(2 * shape2));

    if(x < trunc) break;
    ++i;
  }
  return x;
}

// From https://stackoverflow.com/questions/29724083
uvec arma_setdiff(uvec x, uvec y){

  x = arma::unique(x);
  y = arma::unique(y);

  for (size_t j = 0; j < y.n_elem; j++) {
    uvec q1 = find(x == y[j]);
    if (!q1.empty()) {
      x.shed_row(q1(0));
    }
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
  const uvec& io_idx_input = arma::sort(idx_x);
  std::string message = "C++ indices detected. Unchanged.";
  if (arma::any(io_idx_input - io_idx_cpp)) {
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
  Rcpp::NumericVector set_length_Rcpp, set_length_Rcpp_fact;
  set_length_Rcpp = set_length;
  set_length_Rcpp_fact = Rcpp::factorial(set_length_Rcpp);
  prob /= Rcpp::as<double>(set_length_Rcpp_fact);
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

vec arma_vec_seq(int N) {
  // Creates an arma vector filled with {1,...,N}
  Rcpp::IntegerVector vec_Rcpp = Rcpp::seq(1, N);
  vec v = Rcpp::as<vec>(vec_Rcpp);
  return(v);
}
