#include <RcppArmadillo.h>
#include <cmath>
#include "distances.h"
#include "misc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec propose_augmentation(const vec& ranks, const uvec& indicator){
  vec proposal = ranks;
  proposal(arma::find(indicator == 1)) = arma::shuffle(ranks(arma::find(indicator == 1)));
  return(proposal);
}

void initialize_missing_ranks(mat& rankings, const umat& missing_indicator,
                              const uvec& assessor_missing) {

  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0) {
      continue;
    } else {
      vec rank_vector = rankings.col(i);
      uvec present_inds = arma::find(missing_indicator.col(i) == 0);
      uvec missing_inds = arma::find(missing_indicator.col(i) == 1);
      // Find the available ranks and permute them
      uvec new_ranks = arma::shuffle(arma_setdiff(
        arma::linspace<uvec>(1, rank_vector.size()),
        arma::conv_to<uvec>::from(rank_vector(present_inds))
      ));

      for(unsigned int j = 0; j < missing_inds.size(); ++j){
        rank_vector(missing_inds(j)) = static_cast<double>(arma::as_scalar(new_ranks(j)));
      }
      rankings.col(i) = rank_vector;
    }
  }
}

void update_missing_ranks(mat& rankings, const uvec& current_cluster_assignment,
                          vec& aug_acceptance,
                          const umat& missing_indicator,
                          const uvec& assessor_missing,
                          const vec& alpha, const mat& rho,
                          const std::string& metric){

  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0){
      ++aug_acceptance(i);
      continue;
    }

    // Sample an augmentation proposal
    vec proposal = propose_augmentation(rankings.col(i), missing_indicator.col(i));

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
    }
  }
}
