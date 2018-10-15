#include "RcppArmadillo.h"
#include "misc.h"
#include "distfuns.h"
#include <cmath>


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// Helper to compute the importance sampling smoothed fit
double compute_is_fit(double alpha, arma::vec fit){
  // The partition function
  double logZ = 0;
  int n_items = fit.n_elem;

  for(int i = 0; i < n_items; ++i){
    logZ += std::pow(alpha, i) * fit(i);
  }
  return(logZ);
}




//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param n_items Number of items.
//' @param alpha The value of the alpha parameter.
//' @param cardinalities Number of occurrences for each unique distance.
//' Applicable for footrule and Spearman distance. Defaults to \code{R_NilValue}.
//' @param logz_estimate Precomputed importance sampling fit.
//' @param metric A string. Available options are \code{"footrule"},
//' \code{"kendall"}, \code{"spearman"}, \code{"cayley"}, and \code{"hamming"}.
//' Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @keywords internal
//'
// [[Rcpp::export]]
double get_partition_function(int n_items, double alpha,
                              Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                              Rcpp::Nullable<arma::vec> logz_estimate = R_NilValue,
                              std::string metric = "footrule"){

  if(metric == "footrule") {

    // If cardinalities are defined, we use them. If importance sampling estimates exist,
    // then we use that
    if(cardinalities.isNotNull()){
      arma::vec distances = arma::regspace(0, 2, floor(std::pow(n_items, 2) / 2));
      return std::log(arma::sum(Rcpp::as<arma::vec>(cardinalities) % exp(-alpha * distances / n_items)));
    } else if(logz_estimate.isNotNull()){
      return compute_is_fit(alpha, Rcpp::as<arma::vec>(logz_estimate));
    } else {
      Rcpp::stop("Could not compute partition function");
    }


  } else if (metric == "spearman") {

    if(cardinalities.isNotNull()){
      arma::vec distances = arma::regspace(0, 2, 2 * binomial_coefficient(n_items + 1, 3));
      return std::log(arma::sum(Rcpp::as<arma::vec>(cardinalities) % exp(-alpha * distances / n_items)));
    } else if(logz_estimate.isNotNull()){
      return compute_is_fit(alpha, Rcpp::as<arma::vec>(logz_estimate));
    }


  } else if(metric == "kendall") {

    double res = 0;
    for(int i = 1; i < (n_items + 1); ++i){
      res += std::log( ( 1 - exp(-i * alpha / n_items ) )/(1 - exp(-alpha / n_items )));
    }
    return res;

  } else if(metric == "cayley") {

    double res = 0;

    for(int i = 1; i < n_items; ++i){
      res += std::log( 1 + i * exp(- alpha / n_items ) );
    }
    return res;

  } else if(metric == "hamming"){

    double res = 0;

    for(int i = 0; i < (n_items + 1); ++i){
      res += factorial(n_items) * exp(-alpha) * std::pow((exp(alpha / n_items) - 1), i) / factorial(i);
    }

    return std::log(res);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return 0;
}


//' Asymptotic Approximation of Partition Function
//'
//' Compute the asymptotic approximation of the logarithm of the partition function,
//' using the iteration algorithm of \insertCite{mukherjee2016;textual}{BayesMallows}.
//'
//' @param alpha_vector A numeric vector of alpha values.
//' @param n_items Integer specifying the number of items.
//' @param metric One of \code{"footrule"} and \code{"spearman"}.
//' @param K Integer.
//' @param n_iterations Integer specifying the number of iterations.
//'
//' @return A vector, containing the partition function at each value of alpha.
//' @keywords internal
//'
//' @references \insertAllCited{}
//'
// [[Rcpp::export]]
arma::vec asymptotic_partition_function(arma::vec alpha_vector, int n_items, std::string metric,
                                        int K, int n_iterations){
  // IPFP procedure
  // Initialize a square matrix where each row/column sums to one
  arma::mat A = arma::ones<arma::mat>(K, K) * 1.0 / K;

  // arma::accu sums all elements of the tensor
  // the % operator computes the element-wise product
  double Z0lim = -2 * std::log(K) - arma::accu(A % arma::log(A));

  // Helper matrix
  arma::mat B(K, K);
  for(int i = 0; i < K; ++i){
    for(int j = 0; j < K; ++j){
      if(metric == "footrule"){
        // Because 0 is the first index, this is really (i + 1) - (j + 1) = i - j
        B(i, j) = -std::abs(static_cast<double>(i - j));
      } else if(metric == "spearman"){
        B(i, j) = -std::pow(i - j, 2.0);
      }
    }
  }
  // Divide each element by K
  B = B * 1.0 / K;

  int n_alpha = alpha_vector.n_elem;
  arma::vec result(n_alpha);

  for(int i = 0; i < n_alpha; ++i){
    double alpha = alpha_vector(i);

    A = arma::exp(alpha * B);
    for(int i = 0; i < n_iterations; ++i){
      // Note: We can use 1-norm because the exponential never gets negative
      // Normalize rows
      A = arma::normalise(A, 1, 1);
      // Normalize columns
      A = arma::normalise(A, 1, 0);
    }

    double Zlim = alpha * arma::accu(B % A) - 2 * std::log(K) - arma::accu(A % arma::log(A));

    double Z0 = 0;
    for(int i = 0; i < n_items; ++i){
      Z0 += std::log(static_cast<double>(i + 1));
    }

    result(i) = (Zlim - Z0lim)/K * n_items + Z0;
  }

  return(result);
}




