#ifndef SETDIFF_H
#define SETDIFF_H

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

template <typename arg1, typename arg2>
arma::vec setdiff_template(const arg1& x, const arg2& y, const bool& sort_unique = true){
  Rcpp::NumericVector x_Rcpp, y_Rcpp;
  arma::vec x_y_diff;
  x_Rcpp = x;
  y_Rcpp = y;
  if (sort_unique) {
    x_y_diff = Rcpp::sort_unique(Rcpp::setdiff(x_Rcpp, y_Rcpp));
  } else {
    x_y_diff = Rcpp::setdiff(x_Rcpp, y_Rcpp);
  }

  return x_y_diff;
}
#endif
