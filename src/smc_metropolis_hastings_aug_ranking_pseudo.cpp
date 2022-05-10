#include <RcppArmadillo.h>
#include "smc.h"
#include "misc.h"
#include "setdiff.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Metropolis-Hastings Augmented Ranking (pseudolikelihood)
//' @description Function to perform Metropolis-Hastings for new augmented ranking using the pseudolikelihood augmentation approach
//'
//' @param alpha Numeric value of the scale parameter
//' @param rho Numeric vector specifying the consensus ranking
//' @param n_items Integer is the number of items in a ranking
//' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obtained from calculate_forward_probability function
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return = proposed augmented ranking or current ranking A ranking sequence vector representing proposed augmented ranking for next
//'         iteration of MCMC chain
//' @export
// [[Rcpp::export]]
arma::vec metropolis_hastings_aug_ranking_pseudo(
	const double alpha,
	const arma::vec rho,
	const int n_items,
	const arma::vec partial_ranking,
	const arma::vec current_ranking,
	const std::string metric
) {
  return metropolis_hastings_aug_ranking_both(
    alpha, rho, n_items, partial_ranking, current_ranking, metric, true
  );
}
