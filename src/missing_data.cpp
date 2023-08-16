#include <RcppArmadillo.h>
#include "distances.h"
#include "misc.h"
#include "setdiff.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

uvec propose_augmentation(const uvec& ranks, const uvec& indicator){
  uvec proposal = ranks;
  proposal(find(indicator == 1)) = shuffle(ranks(find(indicator == 1)));
  return proposal;
}

void initialize_missing_ranks(umat& rankings, const umat& missing_indicator,
                              const uvec& assessor_missing) {
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0) {
      continue;
    } else {
      uvec rank_vector = rankings.col(i);
      uvec present_ranks = rank_vector(find(missing_indicator.col(i) == 0));
      uvec missing_inds = find(missing_indicator.col(i) == 1);
      // Find the available ranks and permute them
      vec new_ranks = shuffle(setdiff_template(regspace<vec>(1, rank_vector.n_elem), present_ranks));

      for(unsigned int j = 0; j < missing_inds.n_elem; ++j){
        rank_vector(missing_inds(j)) = new_ranks(j);
      }
      rankings.col(i) = rank_vector;
    }
  }
}

void update_missing_ranks(umat& rankings, const uvec& current_cluster_assignment,
                          vec& aug_acceptance,
                          const umat& missing_indicator,
                          const uvec& assessor_missing,
                          const vec& alpha, const umat& rho,
                          const std::string& metric) {
  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    if(assessor_missing(i) == 0){
      ++aug_acceptance(i);
      continue;
    }

    // Sample an augmentation proposal
    uvec proposal = propose_augmentation(rankings.col(i), missing_indicator.col(i));

    // Draw a uniform random number
    double u = std::log(randu<double>());

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
