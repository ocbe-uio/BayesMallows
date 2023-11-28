#include <RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec normalize_weights(const vec& log_inc_wgt){
  double maxw = max(log_inc_wgt);
  double log_sum_exp = maxw + log(sum(exp(log_inc_wgt - maxw)));
  return exp(log_inc_wgt - log_sum_exp);
}

// Truncated beta distribution
double rtruncbeta(int shape1, int shape2, double trunc = 1) {
  int i = 0;
  double x;
  while(i < 1000){
    x = chi2rnd(2 * shape1);
    x = x / (x + chi2rnd(2 * shape2));

    if(x < trunc) break;
    ++i;
  }
  return x;
}

