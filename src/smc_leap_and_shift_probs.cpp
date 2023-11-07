#include <RcppArmadillo.h>
#include "leapandshift.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List leap_and_shift_probs(const arma::vec rho, const int n_items, const int leap_size = 1) {
  vec rho_proposal{};
  uvec indices{};
  double prob_forward, prob_backward;
  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho, leap_size, false);

  // return(leap_shift_list)
  return Rcpp::List::create(
    Rcpp::Named("rho_prime") = rho_proposal,
    Rcpp::Named("forwards_prob") = prob_forward,
    Rcpp::Named("backwards_prob") = prob_backward
  );
}
