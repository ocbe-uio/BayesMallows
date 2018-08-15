#include <algorithm>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//using namespace Rcpp;

// Function to compute the factorial
// taken from http://www.cplusplus.com/forum/unices/33379/
int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


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


// This function is taken from
// https://stackoverflow.com/questions/29724083/trying-to-write-a-setdiff-function-using-rcpparmadillo-gives-compilation-error
arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {

  std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;

  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));

  return arma::conv_to<arma::uvec>::from(out);
}
