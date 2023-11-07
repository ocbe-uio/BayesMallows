#include <RcppArmadillo.h>
#include "distances.h"
#include "smc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec get_sample_probabilities(
  const vec rho_item_rank,
  const double alpha,
  const vec remaining_set_ranks,
  const int n_items,
  const std::string metric = "footrule"
) {
  // define a set of probs list
  unsigned int num_ranks = remaining_set_ranks.n_elem;
  vec sample_prob_list = zeros(num_ranks);

  // cycle through each item and calculate its specific prob
  for (uword ii = 0; ii < num_ranks; ++ii) {
    const vec item_rank{remaining_set_ranks(ii)};
    const double rank_dist = get_rank_distance(rho_item_rank, item_rank, metric);
    const double sample_prob = -(alpha / n_items) * rank_dist;
    sample_prob_list(ii) = sample_prob;
  }

  return normalize_weights(sample_prob_list);
}
