#include <algorithm>
#include <RcppArmadillo.h>

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
  double u = arma::randu();

  arma::uvec matches = arma::find(arma::cumsum(probs) > u, 1, "first");

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
arma::uvec arma_setdiff(arma::uvec x, arma::uvec y){

  x = arma::unique(x);
  y = arma::unique(y);

  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) {
      x.shed_row(q1(0));
    }
  }
  return x;
}
