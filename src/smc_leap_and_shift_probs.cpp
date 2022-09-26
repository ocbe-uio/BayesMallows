#include <RcppArmadillo.h>
#include "leapandshift.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Leap and Shift Probabilities
//' @description Calculates transition probabilities for proposing a new rho
//' @param rho A ranking sequence
//' @param leap_size Integer specifying the step size of the leap-and-shift
//' proposal distribution.
//' @param n_items Integer is the number of items in a ranking
//' @export
//' @return A list containing:
//' \itemize{
//' \item \code{rho_prime} A ranking sequence proposed consensus ranking
//' \item \code{forwards_prob} Numeric value to account for transition probability from rho to rho_prime
//' \item \code{backwards_prob} Numeric Value to account for the transition probability from \code{rho_prime} to \code{rho}
//' }
//'
//' @keywords internal
//' @examples
//' rho <- c(1, 2, 3, 4, 5, 6)
//' n_items <- 6
//'
//' leap_and_shift_probs(rho, 1, n_items)
//' leap_and_shift_probs(rho, 2, n_items)
//' leap_and_shift_probs(rho, 3, n_items)
//'
// [[Rcpp::export]]
Rcpp::List leap_and_shift_probs(const arma::vec rho, const int leap_size, const int n_items) {
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
