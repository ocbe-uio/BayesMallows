#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export]]
int timesTwo(int x) {
  return x * 2;
}

//' Compute the distance between two rank vectors.
//'
//' @param r1 A vector of ranks.
//' @param r2 A vector of ranks.
//' @param metric A string. Avaiable options are \code{"footrule"},
//' \code{"kendall"}, and \code{"spearman"}. Defaults to \code{"footrule"}.
//' @return A scalar.
//' @export
// [[Rcpp::export]]
double get_rank_distance(const arma::vec r1, const arma::vec r2,
                         const std::string metric = "footrule"){

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
  } else {
    distance = 0;
  }

  return distance;
  // if(metric == "footrule"){
  //   for(int i=0; i<N; i++){
  //     totDist += norm(R.row(i) - rho, 1);
  //   }
  // } else if(metric == "kendall"){
  //   for(int i=0; i<N; i++){
  //     for(int u=0; u<n; u++){
  //       for(int t=0; t<u; t++){
  //         if((R(i,t) > R(i,u)) & (rho(t) < rho(u) ) || ((R(i,t) < R(i,u) ) & (rho(t) > rho(u)))) totDist +=1;
  //       }
  //     }
  //   }
  // } else if(metric == "spearman") {
  //   for(int i=0; i<N; i++){
  //     totDist += pow(norm(R.row(i) - rho, 2),2.0);
  //   }
  // }
  // return totDist;
}
