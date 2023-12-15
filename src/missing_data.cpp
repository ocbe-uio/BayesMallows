#include <RcppArmadillo.h>
#include <algorithm>
#include <Rmath.h>
#include "distances.h"
#include "missing_data.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec setdiff(const vec& x, const vec& y){
  vec xs = sort(x);
  vec ys = sort(y);

  std::vector<double> diff;
  std::set_difference(
    xs.begin(), xs.end(), ys.begin(), ys.end(),
    std::inserter(diff, diff.begin()));

  return conv_to<vec>::from(diff);
}

vec propose_augmentation(const vec& ranks, const uvec& indicator){
  vec proposal = ranks;
  vec mutable_ranks = ranks(find(indicator == 1));
  Rcpp::IntegerVector inds = Rcpp::sample(mutable_ranks.size(),
                                          mutable_ranks.size()) - 1;
  arma::uvec inds_arma = Rcpp::as<uvec>(Rcpp::wrap(inds));
  proposal(find(indicator == 1)) = mutable_ranks(inds_arma);
  return proposal;
}

void set_up_missing(arma::mat& rankings, arma::umat& missing_indicator) {
  rankings.replace(datum::nan, 0);
  missing_indicator = conv_to<umat>::from(rankings);
  missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
}

void initialize_missing_ranks(mat& rankings, const umat& missing_indicator) {
  int n_assessors = rankings.n_cols;

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
    rankings.col(i) = rank_vector;
  }
}

vec make_new_augmentation(const vec& rankings, const uvec& missing_indicator,
                          const double& alpha, const vec& rho,
                          const std::string& metric, bool pseudo) {
  double log_hastings_correction = 0;
  vec proposal{};
  if(pseudo) {
    uvec a = find(missing_indicator == 1);
    Rcpp::IntegerVector b = Rcpp::sample(a.size(), a.size()) - 1;
    uvec unranked_items = a(Rcpp::as<uvec>(Rcpp::wrap(b)));

    auto pprop = make_pseudo_proposal(
      unranked_items, rankings, alpha, rho, metric, true
    );

    auto bprop = make_pseudo_proposal(
      unranked_items, rankings, alpha, rho, metric, false);

    proposal = pprop.rankings;
    log_hastings_correction = -std::log(pprop.probability) +
      std::log(bprop.probability);
  } else {
    proposal = propose_augmentation(rankings, missing_indicator);
  }

  double u = std::log(R::runif(0, 1));
  int n_items = rho.n_elem;

  double ratio = -alpha / n_items *
    (get_rank_distance(proposal, rho, metric) -
    get_rank_distance(rankings, rho, metric)) + log_hastings_correction;

  if(ratio > u){
    return proposal;
  } else {
    return rankings;
  }
}

PseudoProposal make_pseudo_proposal(
    uvec unranked_items, vec rankings, const double& alpha, const vec& rho,
    const std::string metric, const bool forward
) {
  int n_items = rankings.n_elem;
  double prob = 1;
  while(unranked_items.n_elem > 0) {
    vec available_rankings = rankings(unranked_items);
    int item_to_rank = unranked_items(0);

    vec log_numerator(available_rankings.n_elem);
    for(size_t ll{}; ll < available_rankings.n_elem; ll++) {
      log_numerator(ll) = -alpha / n_items *
        get_rank_distance(rho(span(item_to_rank)), available_rankings(span(ll)), metric);
    }
    vec sample_probs = normalise(exp(log_numerator), 1);
    if(forward) {
      ivec ans(sample_probs.size());
      R::rmultinom(1, sample_probs.begin(), sample_probs.size(), ans.begin());
      rankings(span(item_to_rank)) = available_rankings(find(ans == 1));
    }

    int ranking_chosen = as_scalar(find(rankings(item_to_rank) == available_rankings));

    prob *= sample_probs(ranking_chosen);
    if(available_rankings.n_elem <= 1) break;
    unranked_items = unranked_items.subvec(1, available_rankings.n_elem - 1);

    rankings(unranked_items) = setdiff(
      available_rankings, available_rankings(span(ranking_chosen)));
  }

  return PseudoProposal(rankings, prob);
}
