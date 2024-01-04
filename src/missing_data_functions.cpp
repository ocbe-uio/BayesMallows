#include <RcppArmadillo.h>
#include <algorithm>
#include "distances.h"
#include "missing_data.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec setdiff(const vec& x, const vec& y) noexcept {
  vec xs = sort(x);
  vec ys = sort(y);

  std::vector<double> diff;
  std::set_difference(
    xs.begin(), xs.end(), ys.begin(), ys.end(),
    std::inserter(diff, diff.begin()));

  return conv_to<vec>::from(diff);
}

arma::umat set_up_missing(const Data& dat) noexcept {
  if(!dat.any_missing) return arma::umat{};
  arma::umat missing_indicator = conv_to<umat>::from(dat.rankings);
  missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
  return missing_indicator;
}

mat initialize_missing_ranks(const mat& rankings, const umat& missing_indicator) {
  int n_assessors = rankings.n_cols;
  mat init_rank = rankings;

  for(int i = 0; i < n_assessors; ++i){
    vec rank_vector = rankings.col(i);
    vec present_ranks = rank_vector(find(missing_indicator.col(i) == 0));
    uvec missing_inds = find(missing_indicator.col(i) == 1);
    vec a = setdiff(
      regspace<vec>(1, rank_vector.n_elem), present_ranks);

    Rcpp::IntegerVector inds = Rcpp::sample(a.size(), a.size()) - 1;
    uvec inds_arma = Rcpp::as<uvec>(Rcpp::wrap(inds));
    vec new_ranks = a.elem(inds_arma);
    rank_vector(missing_inds) = new_ranks;
    init_rank.col(i) = rank_vector;
  }
  return init_rank;
}

vec make_new_augmentation(const vec& rankings, const uvec& missing_indicator,
                          const double& alpha, const vec& rho,
                          const std::unique_ptr<Distance>& distfun,
                          const std::unique_ptr<Distance>& pseudo_aug_distance,
                          double& log_aug_prob) {
  double log_hastings_correction = 0;
  RankProposal pprop{};
  if(pseudo_aug_distance == nullptr) {
    pprop = make_uniform_proposal(rankings, missing_indicator);
  } else {
    pprop = make_pseudo_proposal(
      rankings, missing_indicator, alpha, rho, pseudo_aug_distance);
    log_hastings_correction = -std::log(pprop.probability) + log_aug_prob;
  }

  double u = std::log(R::runif(0, 1));
  int n_items = rho.n_elem;

  double newdist = distfun->d(pprop.rankings, rho, pprop.mutated_items);
  double olddist = distfun->d(rankings, rho, pprop.mutated_items);
  double ratio = -alpha / n_items * (newdist - olddist) +
    log_hastings_correction;

  if(ratio > u){
    log_aug_prob = std::log(pprop.probability);
    return pprop.rankings;
  } else {
    return rankings;
  }
}

RankProposal make_uniform_proposal(const vec& ranks, const uvec& indicator) noexcept {
  vec proposal = ranks;
  uvec missing_inds = find(indicator == 1);
  vec mutable_ranks = ranks(missing_inds);
  Rcpp::IntegerVector inds = Rcpp::sample(mutable_ranks.size(),
                                          mutable_ranks.size()) - 1;
  uvec inds_arma = Rcpp::as<uvec>(Rcpp::wrap(inds));
  proposal(missing_inds) = mutable_ranks(inds_arma);
  return RankProposal(proposal, 1, missing_inds);
}

RankProposal make_pseudo_proposal(
    vec ranks, const uvec& indicator, const double& alpha, const vec& rho,
    const std::unique_ptr<Distance>& distfun
) noexcept {
  uvec missing_inds = find(indicator == 1);
  vec a = Rcpp::runif(missing_inds.size(), 0, 1);
  uvec b = sort_index(a);
  uvec unranked_items = missing_inds(b);

  int n_items = ranks.n_elem;
  double prob = 1;
  while(unranked_items.n_elem > 0) {
    vec available_rankings = ranks(unranked_items);
    int item_to_rank = unranked_items(0);

    double rho_for_item = rho(item_to_rank);
    vec log_numerator = -alpha / n_items *
      distfun->d(available_rankings, rho_for_item);

    vec sample_probs = normalise(exp(log_numerator), 1);
    ivec ans(sample_probs.size());
    R::rmultinom(1, sample_probs.begin(), sample_probs.size(), ans.begin());
    ranks(span(item_to_rank)) = available_rankings(find(ans == 1));

    int ranking_chosen = as_scalar(find(ranks(item_to_rank) == available_rankings));

    prob *= sample_probs(ranking_chosen);
    if(available_rankings.n_elem <= 1) break;
    unranked_items = unranked_items.subvec(1, available_rankings.n_elem - 1);

    ranks(unranked_items) = setdiff(
      available_rankings, available_rankings(span(ranking_chosen)));
  }

  return RankProposal(ranks, prob, missing_inds);
}
