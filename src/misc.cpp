#include <RcppArmadillo.h>
#include <algorithm>

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


// Temporary development
// [[Rcpp::export]]
arma::mat missing_data(arma::mat R){
  int n_items = R.n_rows;
  int n_assessors = R.n_cols;

  // Declare indicator matrix of missing ranks, and fill it with zeros
  arma::mat missing_indicator = arma::zeros<arma::mat>(n_items, n_assessors);

  // Number of missing items per assessor
  arma::vec assessor_missing = arma::zeros<arma::vec>(n_assessors);

  for(int i = 0; i < n_assessors; ++i){

    // Update the missingness indicator
    for(int j = 0; j < n_items; ++j){
      if(!arma::is_finite(R(j, i))){
        // Set the indicator matrix element to 1
        missing_indicator(j, i) = 1;

        // Increment the missingness count for the assessor
        ++assessor_missing(i);
      }
    }

    // Go to next step in loop if assessor does not have missing values
    if(assessor_missing(i) == 0) continue;

    // Defining a number of helper variables
    arma::uvec ranked_inds = arma::find(missing_indicator.col(i) == 0);
    arma::uvec taken_ranks = arma::conv_to<arma::uvec>::from(R.col(i)).elem(ranked_inds);
    arma::uvec all_ranks = arma::regspace<arma::uvec>(1, 1, n_items);
    arma::uvec all_inds = arma::regspace<arma::uvec>(0, 1, n_items - 1);

    // Find the available ranks
    arma::uvec available_ranks = std_setdiff(all_ranks, taken_ranks);
    arma::uvec available_inds = std_setdiff(all_inds, ranked_inds);

    // Fill in available ranks
    arma::vec new_ranks = R.col(i);
    // Adding randomness with shuffle
    new_ranks(available_inds) = arma::shuffle(arma::conv_to<arma::vec>::from(available_ranks));
    R.col(i) = new_ranks;



    // miss = find(R.row(i) < 0); // which elements to augment
    // tmp = trans(R.row(i));
    // tmp(miss) = shuffle(linspace(Rmiss(i), n, n-Rmiss(i)+1));
    // Raug.row(i) = trans(tmp);

  }

  return R;
}
