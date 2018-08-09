#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

void define_missingness(arma::mat& missing_indicator, arma::vec& assessor_missing,
                        const arma::mat& R,
                        const int& n_items, const int& n_assessors){
  for(int i = 0; i < n_assessors; ++i){
    for(int j = 0; j < n_items; ++j){
      if(!arma::is_finite(R(j, i))){
        missing_indicator(j, i) = 1;
        ++assessor_missing(i);
      }
    }
  }

}
