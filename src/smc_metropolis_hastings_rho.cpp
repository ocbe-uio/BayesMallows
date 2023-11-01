#include <RcppArmadillo.h>
#include "smc.h"
#include "leapandshift.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec metropolis_hastings_rho(
	const double alpha,
	const int n_items,
	const arma::mat rankings,
	const arma::vec rho,
	const std::string metric = "footnote",
	const int leap_size = 1
) {
  // create new potential consensus ranking
  vec rho_proposal{};
  uvec indices{};
  double prob_forward, prob_backward;
  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho, leap_size, false);

  // evaluate the log-likelihood with current rankings
  const double mallows_loglik_curr = get_exponent_sum(alpha, rho, n_items, rankings, metric);
  const double mallows_loglik_prop = get_exponent_sum(alpha, rho_proposal, n_items, rankings, metric);

  // calculate acceptance probability
  const double& loga = std::log(prob_backward) - std::log(prob_forward) +
    mallows_loglik_prop - mallows_loglik_curr;

  // determine whether to accept or reject proposed rho and return now consensus ranking
  const double& p = R::runif(0, 1);
  if (std::log(p) <= loga) {
    return(rho_proposal);
  } else {
    return(rho);
  }
}
