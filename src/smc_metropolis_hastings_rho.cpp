#include <RcppArmadillo.h>
#include "smc.h"
#include "leapandshift.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Rho
//' @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!
//' @inheritParams get_exponent_sum
//' @param leap_size Integer specifying the step size of the leap-and-shift
//' proposal distribution.
//' @export
//' @author Anja Stein
//' @examples
//' rho <- t(c(1,2,3,4,5,6))
//' alpha <- 2
//' metric <- "footrule"
//' n_items <- 6
//'
//' metropolis_hastings_rho(
//' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
//' 	rho = rho, leap_size = 1
//' )
//'
//' metropolis_hastings_rho(
//' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
//' 	rho = rho, leap_size = 2
//' )
//'
//' metropolis_hastings_rho(
//' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
//' 	rho = rho, leap_size = 3
//' )
//'
//' rankings <- sample_mallows(
//'  rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
//' )
//' metropolis_hastings_rho(
//' 	alpha = alpha, n_items = n_items, rankings = rankings, metric = metric,
//' 	rho = rho, leap_size = 1
//' )
//'
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
