#ifndef SETDIFF_H
#define SETDIFF_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

template <typename T1, typename T2>
arma::vec setdiff(const T1& x, const T2& y){
  Rcpp::NumericVector x_Rcpp, y_Rcpp;
  arma::vec x_y_diff;
  x_Rcpp = x;
  y_Rcpp = y;
  x_y_diff = Rcpp::sort_unique(Rcpp::setdiff(x_Rcpp, y_Rcpp));

  return x_y_diff;
}
#endif
