#include <RcppArmadillo.h>
#include <algorithm>
#include "distances.h"
#include "missing_data.h"
#include "setdiff.h"

using namespace arma;

mat initialize_missing_ranks(mat rankings, const umat& missing_indicator) {
  int n_assessors = rankings.n_cols;

  for(int i = 0; i < n_assessors; ++i){
    vec rank_vector = rankings.col(i);
    vec present_ranks = rank_vector(find(missing_indicator.col(i) == 0));
    uvec missing_inds = find(missing_indicator.col(i) == 1);
    vec a = setdiff(
      regspace<vec>(1, rank_vector.n_elem), present_ranks);

    ivec inds = Rcpp::sample(a.size(), a.size()) - 1;
    vec new_ranks = a.elem(conv_to<uvec>::from(inds));
    rank_vector(missing_inds) = new_ranks;
    rankings.col(i) = rank_vector;
  }
  return rankings;
}


vec make_new_augmentation_common(const vec& rankings, double alpha, const vec& rho,
                                 const std::unique_ptr<Distance>& distfun,
                                 const RankProposal& pprop, double& log_aug_prob,
                                 const std::string& error_model = "none", double theta = 0.0) {

  double log_hastings_correction = -std::log(pprop.prob_forward) + log_aug_prob;

  double newdist = distfun->d(pprop.rankings, rho, pprop.mutated_items);
  double olddist = distfun->d(rankings, rho, pprop.mutated_items);
  double ratio = -alpha / rho.n_elem * (newdist - olddist) + log_hastings_correction;

  if(error_model != "none") {
    ratio += pprop.g_diff * std::log(theta / (1 - theta));
  }

  if(ratio > std::log(R::runif(0, 1))){
    log_aug_prob = std::log(pprop.prob_forward);
    return pprop.rankings;
  } else {
    return rankings;
  }
}

vec make_new_augmentation(const vec& rankings, const uvec& missing_indicator,
                          double alpha, const vec& rho,
                          const std::unique_ptr<Distance>& distfun,
                          const std::unique_ptr<PartialProposal>& partial_aug_prop,
                          double& log_aug_prob) {

  RankProposal pprop = partial_aug_prop->propose(rankings, missing_indicator, alpha, rho);
  return make_new_augmentation_common(rankings, alpha, rho, distfun, pprop, log_aug_prob);
}

vec make_new_augmentation(const vec& rankings, double alpha, const vec& rho,
                          double theta,
                          const std::unique_ptr<Distance>& distfun,
                          const std::unique_ptr<PairwiseProposal>& pairwise_aug_prop,
                          double& log_aug_prob, const doubly_nested& items_above,
                          const doubly_nested& items_below, const std::string& error_model) {

  RankProposal pprop = pairwise_aug_prop->propose(rankings, items_above, items_below);

  return make_new_augmentation_common(rankings, alpha, rho, distfun, pprop, log_aug_prob, error_model, theta);
}
