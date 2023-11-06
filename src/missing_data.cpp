#include <RcppArmadillo.h>
#include "distances.h"
#include "misc.h"
#include "setdiff.h"
#include "missing_data.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec propose_augmentation(const vec& ranks, const uvec& indicator){
  vec proposal = ranks;
  proposal(find(indicator == 1)) = shuffle(ranks(find(indicator == 1)));
  return proposal;
}

void initialize_missing_ranks(mat& rankings, const umat& missing_indicator) {
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    vec rank_vector = rankings.col(i);
    vec present_ranks = rank_vector(find(missing_indicator.col(i) == 0));
    uvec missing_inds = find(missing_indicator.col(i) == 1);
    // Find the available ranks and permute them
    vec new_ranks = shuffle(setdiff_template(regspace<vec>(1, rank_vector.n_elem), present_ranks));

    for(unsigned int j = 0; j < missing_inds.n_elem; ++j){
      rank_vector(missing_inds(j)) = new_ranks(j);
    }
    rankings.col(i) = rank_vector;

  }
}

void update_missing_ranks(mat& rankings, const uvec& current_cluster_assignment,
                          const umat& missing_indicator,
                          const vec& alpha, const mat& rho,
                          const std::string& metric) {
  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){

    int cluster = current_cluster_assignment(i);

    rankings.col(i) = make_new_augmentation(
      rankings.col(i), missing_indicator.col(i),
      alpha(cluster), rho.col(cluster), metric
    );

  }
}

vec make_new_augmentation(const vec& rankings, const uvec& missing_indicator,
                          const double& alpha, const vec& rho,
                          const std::string& metric) {
  // Sample an augmentation proposal
  vec proposal = propose_augmentation(rankings, missing_indicator);

  // Draw a uniform random number
  double u = std::log(randu<double>());

  int n_items = rho.n_elem;

  double ratio = -alpha / n_items *
    (get_rank_distance(proposal, rho, metric) -
    get_rank_distance(rankings, rho, metric));

  if(ratio > u){
    return proposal;
  } else {
    return rankings;
  }
}
