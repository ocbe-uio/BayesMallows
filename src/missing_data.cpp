#include <RcppArmadillo.h>
#include <cmath>
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec propose_augmentation(const arma::vec& ranks, const arma::vec& indicator){
  arma::vec proposal = ranks;
  proposal(arma::find(indicator == 1)) = arma::shuffle(ranks(arma::find(indicator == 1)));
  return(proposal);
}

void initialize_missing_ranks(arma::mat& rankings, const arma::mat& missing_indicator,
                              const arma::vec& assessor_missing) {

  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0) {
      continue;
    } else {
      rankings.col(i) = propose_augmentation(rankings.col(i), missing_indicator.col(i));
    }
  }
}

void update_missing_ranks(arma::mat& rankings, const arma::uvec& current_cluster_assignment,
                          arma::vec& aug_acceptance,
                          const arma::mat& missing_indicator,
                          const arma::vec& assessor_missing,
                          const arma::vec& alpha, const arma::mat& rho,
                          const std::string& metric,
                          bool& augmentation_accepted){

  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0){
      ++aug_acceptance(i);
      continue;
    }

    // Sample an augmentation proposal
    arma::vec proposal = propose_augmentation(rankings.col(i), missing_indicator.col(i));

    // Draw a uniform random number
    double u = std::log(arma::randu<double>());

    // Find which cluster the assessor belongs to
    int cluster = current_cluster_assignment(i);

    double ratio = -alpha(cluster) / n_items *
      (get_rank_distance(proposal, rho.col(cluster), metric) -
      get_rank_distance(rankings.col(i), rho.col(cluster), metric));

    if(ratio > u){
      rankings.col(i) = proposal;
      ++aug_acceptance(i);
      augmentation_accepted = true;
    } else {
      augmentation_accepted = false;
    }
  }
}
