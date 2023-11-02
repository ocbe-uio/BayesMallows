#include <RcppArmadillo.h>
#include "smc.h"
#include "leapandshift.h"
#include "distances.h"

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
  int n_assessors = rankings.n_cols;
  vec rho_proposal{};
  uvec indices{};
  double prob_forward, prob_backward;
  leap_and_shift(rho_proposal, indices, prob_backward, prob_forward,
                 rho, leap_size, false);

  // evaluate the log-likelihood with current rankings
  double dist_new = rank_dist_sum(rankings.rows(indices), rho_proposal(indices), metric, ones(n_assessors));
  double dist_old = rank_dist_sum(rankings.rows(indices), rho(indices), metric, ones(n_assessors));

  // calculate acceptance probability
  double ratio = - alpha / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  double u = std::log(randu<double>());

  if(ratio > u){
    return(rho_proposal);
  } else {
    return(rho);
  }

}
