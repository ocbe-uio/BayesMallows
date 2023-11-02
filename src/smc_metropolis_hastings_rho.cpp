#include <RcppArmadillo.h>
#include "smc.h"
#include "leapandshift.h"
#include "distances.h"
#include "parameterupdates.h"

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

  int n_assessors = rankings.n_cols;
  return(make_new_rho(rho, rankings, alpha, leap_size, metric, ones(n_assessors)));

}
