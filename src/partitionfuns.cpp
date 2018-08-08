#include "RcppArmadillo.h"
#include "misc.h"
#include "distfuns.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// Helper to compute the importance sampling smoothed fit
double compute_is_fit(double alpha, arma::vec fit){
  // The partition function
  double logZ = 0;
  int n = fit.n_elem;

  for(int i = 0; i < n; ++i){
    logZ += pow(alpha, i) * fit(i);
  }
  return(logZ);
}




//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param n Number of items.
//' @param alpha The value of the alpha parameter.
//' @param cardinalities Number of occurences for each unique distance.
//' Applicable for footrule and Spearman distance. Defaults to \code{R_NilValue}.
//' @param is_fit Precomputed importance sampling fit.
//' @param metric A string. Available options are \code{"footrule"},
//' \code{"kendall"}, \code{"spearman"}, \code{"cayley"}, and \code{"hamming"}.
//' Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
// [[Rcpp::export]]
double get_partition_function(int n, double alpha,
                              Rcpp::Nullable<arma::vec> cardinalities = R_NilValue,
                              Rcpp::Nullable<arma::vec> is_fit = R_NilValue,
                              std::string metric = "footrule"){

  if(metric == "footrule") {

    // If cardinalities are defined, we use them. If importance sampling estimates exist,
    // then we use that
    if(cardinalities.isNotNull()){
      arma::vec distances = arma::regspace(0, 2, floor(pow(n, 2) / 2));
      return log(arma::sum(Rcpp::as<arma::vec>(cardinalities) % exp(-alpha * distances / n)));
    } else if(is_fit.isNotNull()){
      return compute_is_fit(alpha, Rcpp::as<arma::vec>(is_fit));
    } else {
      Rcpp::stop("Could not compute partition function");
    }


  } else if (metric == "spearman") {

    if(cardinalities.isNotNull()){
      arma::vec distances = arma::regspace(0, 2, 2 * binomial_coefficient(n + 1, 3));
      return log(arma::sum(Rcpp::as<arma::vec>(cardinalities) % exp(-alpha * distances / n)));
    } else if(is_fit.isNotNull()){
      return compute_is_fit(alpha, Rcpp::as<arma::vec>(is_fit));
    }


  } else if(metric == "kendall") {

    double res = 0;
    for(int i = 1; i < (n + 1); ++i){
      res += log( ( 1 - exp(-i * alpha / n ) )/(1 - exp(-alpha / n )));
    }
    return res;

  } else if(metric == "cayley") {

    double res = 0;

    for(int i = 1; i < n; ++i){
      res += log( 1 + i * exp(- alpha / n ) );
    }
    return res;

  } else if(metric == "hamming"){

    double res = 0;

    for(int i = 0; i < (n + 1); ++i){
      res += factorial(n) * exp(-alpha) * pow((exp(alpha / n) - 1), i) / factorial(i);
    }

    return log(res);

  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return 0;
}


