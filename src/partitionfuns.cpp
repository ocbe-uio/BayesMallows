#include "RcppArmadillo.h"
#include "misc.h"
#include "distfuns.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


//' Compute the logarithm of the partition function for a Mallows rank model.
//'
//' @param n Number of items.
//' @param alpha The value of the alpha parameter.
//' @param cardinalities Number of occurences for each unique distance.
//' Applicable for footrule and Spearman distance.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, \code{"spearman"}, \code{"cayley"}, and \code{"hamming"}.
//' Defaults to \code{"footrule"}.
//' @return A scalar, the logarithm of the partition function.
//' @export
// [[Rcpp::export]]
double get_partition_function(int n, double alpha, arma::vec cardinalities,
                              std::string metric = "footrule"){

  if(metric == "footrule") {

    arma::vec distances = arma::regspace(0, 2, floor(pow(n, 2) / 2));
    return log(arma::sum(cardinalities % exp(-alpha * distances / n)));

  } else if (metric == "spearman") {
    arma::vec distances = arma::regspace(0, 2, 2 * binomial_coefficient(n + 1, 3));
    return log(arma::sum(cardinalities % exp(-alpha * distances / n)));
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

}

