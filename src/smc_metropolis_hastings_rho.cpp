#include "RcppArmadillo.h"
#include "smc.h"
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Rho
//' @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!
//' @inheritParams get_mallows_loglik
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
	double alpha,
	int n_items,
	arma::mat rankings,
	std::string metric,
	arma::vec rho,
	int leap_size
) {
  // create new potential consensus ranking
  Rcpp::List kernel = leap_and_shift_probs(rho, leap_size, n_items);

  // output from leap-and-shift is of the following
  arma::vec rho_prime = Rcpp::as<arma::vec>(kernel["rho_prime"]);
  double forwards_prob = Rcpp::as<double>(kernel["forwards_prob"]); // rho_prime|rho
  double backwards_prob = Rcpp::as<double>(kernel["backwards_prob"]); // rho|rho_prime

  // evaluate the log-likelihood with current rankings
  double mallows_loglik_curr = get_mallows_loglik(alpha, rho, n_items, rankings, metric);
  double mallows_loglik_prop = get_mallows_loglik(alpha, rho_prime, n_items, rankings, metric);

  // calculate acceptance probability
  double loga = std::log(backwards_prob) - std::log(forwards_prob) + mallows_loglik_prop - mallows_loglik_curr;

  // determine whether to accept or reject proposed rho and return now consensus ranking
  double p = Rcpp::as<double>(Rcpp::runif(1, 0, 1));
  if (std::log(p) <= loga) {
    return(rho_prime);
  } else {
    return(rho);
  }
}
