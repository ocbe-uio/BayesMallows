#include "RcppArmadillo.h"
#include "math.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

struct summation_sequence {
  arma::vec distances;
  arma::vec cardinalities;
};


//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param n The number of items, a positive integer.
//' @param alpha The value of the alpha parameter.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
// [[Rcpp::export]]
double get_partition_functions(int n, double alpha, std::string metric = "footrule"){

  if (n <= 0){
    Rcpp::stop("n must be a positive integer");
  }

  double log_z_n = 0;

  if(metric == "footrule") {



  } else {
    Rcpp::stop("Inadmissible value of metric.");
  }

  return log_z_n;
}


summation_sequence get_summation_sequences(int n, std::string metric = "footrule") {

  summation_sequence result;

  if((metric == "footrule") & (n < 51)){
      result.distances = arma::linspace(0, floor(pow(n, 2) / 2.0));
    // must figure out how to get the cardinalities here. The current linspace is just a dummy.
      result.cardinalities = arma::linspace(0, 10);
  }

  return result;
}

