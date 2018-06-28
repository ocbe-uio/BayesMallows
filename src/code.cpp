#include "RcppArmadillo.h"
#include "math.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// We put all the function declarations here for clarity
int binomial_coefficient(int n, int k);
double get_rank_distance(arma::vec, arma::vec, std::string);
double get_partition_function(double, Rcpp::List, std::string);
Rcpp::List run_mcmc(arma::Mat<int>, arma::vec, std::string);
arma::vec get_summation_distances(int, arma::vec, std::string);

//' Worker function for computing the posterior distribtuion.
//'
//' @param R A set of complete rankings.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function.
//' TODO: Make this argument optional.
//' @param metric The distance metric to use.
//' @param nmc Number of Monte Carlo samples.
//' @param L Leap-and-shift step size.
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::Mat<int> R, arma::vec cardinalities,
                    std::string metric = "footrule", int nmc = 10,
                    int L = 1){

  // The number of items ranked
  int n = R.n_cols;

  // First we find the vector of distances used to compute the partition function
  arma::vec distances = get_summation_distances(n, cardinalities, metric);

  // Declare the matrix to hold the latent ranks
  arma::mat rho = arma::zeros<arma::mat>(nmc, n);

  // Set the initial latent rank value
  rho.row(0) = arma::linspace<arma::rowvec>(1, n, n);

  // Declare the vector to hold the scaling parameter alpha
  arma::vec alpha = arma::zeros<arma::vec>(nmc);

  // Set the initial alpha value
  alpha(0) = 1;

  for(int t = 1; t < nmc; ++t){
    rho.row(t) = rho.row(t - 1);
    alpha(t) = alpha(t - 1);
  }

  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("alpha") = alpha
  );
}

//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param alpha The value of the alpha parameter.
//' @param summation_sequences List with elements \code{distances} and
//' \code{cardinalities}, both of type arma::vec.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
// [[Rcpp::export]]
double get_partition_function(double alpha, Rcpp::List summation_sequences,
                              std::string metric = "footrule"){



  double log_z_n = 0;

  if(metric == "footrule") {
    log_z_n = 1;


  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return log_z_n;
}

//' Get the distances for computing the partition function given
//' the cardinalities.
//'
//' @param alpha The value of the alpha parameter.
//' @param summation_sequences List with elements \code{distances} and
//' \code{cardinalities}, both of type arma::vec.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
// [[Rcpp::export]]
arma::vec get_summation_distances(int n, arma::vec cardinalities,
                                   std::string metric = "footrule") {

  arma::vec distances;

  if(metric == "footrule"){

    if(n > 50) Rcpp::stop("n > 50 currently not supported for footrule");

    int max = floor(pow(n, 2) / 2.0);
    // Sequence from 0 to max with increment 1
    distances = arma::linspace(0, max, max + 1);

  } else if(metric == "spearman") {

    if(n > 13) Rcpp::stop("n > 13 currently not supported for Spearman distance");

    int max = 2 * binomial_coefficient(n, 3);
    distances = arma::linspace(0, max, max + 1);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return distances;
}


//' Compute the distance between two rank vectors.
//'
//' @param r1 A vector of ranks.
//' @param r2 A vector of ranks.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar.
//' @details Note that the Spearman distance is the squared L2 norm, whereas
//' the footrule distance is the L1 norm.
//' @export
// [[Rcpp::export]]
double get_rank_distance(arma::vec r1, arma::vec r2, std::string metric = "footrule"){

  if (r1.n_elem != r2.n_elem){
    Rcpp::Rcout << "r1 and r2 must have the same length" << std::endl;
    exit(1);
  }
  int n = r1.n_elem;
  double distance = 0;

  if(metric == "footrule") {

    // Footrule is the one-norm
    distance = arma::norm(r1 - r2, 1);

  } else if(metric == "kendall") {

    // Need loops to compute Kendall distance

    for(int i = 0; i < n; ++i){
      for(int j = 0; j < i; ++j){
        if(((r1(j) > r1(i)) & (r2(j) < r2(i)) ) || ((r1(j) < r1(i)) & (r2(j) > r2(i)))) {
          distance += 1;
        }
      }
    }

  } else if (metric == "spearman") {

    // Spearman distance is the squared L2 norm
    distance = pow(arma::norm(r1 - r2, 2), 2);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return distance;
}


int binomial_coefficient(int n, int k){
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
