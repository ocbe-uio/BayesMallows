#include <RcppArmadillo.h>
#include <cmath>

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
//' @examples
//' rho <- c(1, 2, 3, 4, 5, 6)
//' n_items <- 6
//'
//' leap_and_shift_probs(rho, 1, n_items)
//' leap_and_shift_probs(rho, 2, n_items)
//' leap_and_shift_probs(rho, 3, n_items)
//' @author Anja Stein
//'
// [[Rcpp::export]]
Rcpp::List leap_and_shift_probs(const arma::vec rho, const int leap_size, const int n_items) {

  // draw one u uniformly from {1,...,n} to use as index for rho
  int u = randi<int>(distr_param(0, n_items - 1));

  // define set of integers S, the support set for sampling new ranks
  const int rho_minus_leap = rho(u) - leap_size;
  const int rho_plus_leap = rho(u) + leap_size;
  const int low_bd = std::max(1, rho_minus_leap);
  const int max_bd = std::min(n_items, rho_plus_leap);
  ivec S = regspace<ivec>(low_bd, max_bd);
  S = S.elem(find(S != rho(u)));

  // draw a random index ind from the elements of S
  int ind = randi<int>(distr_param(0, S.n_elem - 1));
  // store the value r in S corresponding to the index ind
  int r = S(ind);

  // Create leap step
  vec rho_star = rho;
  rho_star(u) = r; // replace u-th entry with r

  // here, two elements are the same so we need to shift element and replace the repeated r with u
  const int& delta = rho_star(u) - rho(u);
  ivec rho_prime = zeros<ivec>(n_items);

  // Define support set for ranks rho_star[u] can leap to
  const int rho_star_minus_leap = rho_star(u) - leap_size;
  const int rho_star_plus_leap = rho_star(u) + leap_size;
  const int S_star_min = std::max(1, rho_star_minus_leap);
  const int S_star_max = std::min(n_items, rho_star_plus_leap);
  ivec S_star;
  S_star = Rcpp::seq(S_star_min, S_star_max);
  S_star = S_star.elem(find(S_star != rho_star(u)));

  // calculate forward and backwards probabilities
  vec forwards_prob, backwards_prob;
  if (std::abs(rho_star(u) - rho(u)) == 1) {
    // p(proposed|current)
    forwards_prob = 1.0 / (n_items * S.n_elem) + 1.0 / (n_items * S_star.n_elem);
    // p(current|proposed)
    backwards_prob = forwards_prob;
  } else {
    // p(proposed|current)
    forwards_prob = 1.0 / (n_items * S.n_elem);
    // p(current|proposed)
    backwards_prob = 1.0 / (n_items * S_star.n_elem);
  }

  // shift step
  for (int i = 0; i < n_items; ++i) {
    if (rho(i) == rho(u)) {
      rho_prime(i) = rho_star(u);
    } else if ((rho(u) < rho(i)) && (rho(i) <= rho_star(u)) && (delta > 0)) {
      rho_prime(i) = rho(i) - 1;
    } else if ((rho(u) > rho(i)) && (rho(i) >= rho_star(u)) && (delta < 0)) {
      rho_prime(i) = rho(i) + 1;
    } else {
      rho_prime(i) = rho(i);
    }
  }

  // return(leap_shift_list)
  return Rcpp::List::create(
    Rcpp::Named("rho_prime") = rho_prime,
    Rcpp::Named("forwards_prob") = forwards_prob,
    Rcpp::Named("backwards_prob") = backwards_prob
  );
}
