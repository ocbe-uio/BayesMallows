#include <RcppArmadillo.h>
#include "setdiff.h"
#include "distances.h"
#include "missing_data.h"
#include "misc.h"
#include "sample.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

vec propose_augmentation(const vec& ranks, const uvec& indicator){
  vec proposal = ranks;
  proposal(find(indicator == 1)) = shuffle(ranks(find(indicator == 1)));
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
    // Find the available ranks and permute them
    vec new_ranks = shuffle(setdiff(regspace<vec>(1, rank_vector.n_elem), present_ranks));

    rank_vector(missing_inds) = new_ranks;
    rankings.col(i) = rank_vector;
  }
}

void update_missing_ranks(mat& rankings, const uvec& current_cluster_assignment,
                          const umat& missing_indicator,
                          const vec& alpha, const mat& rho,
                          const std::string& metric) {

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
                          const std::string& metric, bool pseudo) {
  double log_hastings_correction = 0;
  vec proposal{};
  // Sample an augmentation proposal
  if(pseudo) {
    uvec unranked_items = shuffle(find(missing_indicator == 1));
    Rcpp::List pprop = make_pseudo_proposal(
      unranked_items, rankings, alpha, rho, metric, true
    );

    Rcpp::List bprop = make_pseudo_proposal(
      unranked_items, rankings, alpha, rho, metric, false);
    double bprob = bprop["probability"];

    vec ar = pprop["proposal"];
    proposal = ar;
    double prob = pprop["probability"];

    log_hastings_correction = -std::log(prob) + std::log(bprob);

  } else {
    proposal = propose_augmentation(rankings, missing_indicator);
  }


  // Draw a uniform random number
  double u = std::log(randu<double>());

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


Rcpp::List make_pseudo_proposal(
    uvec unranked_items, vec rankings, const double& alpha, const vec& rho,
    const std::string metric, const bool forward
) {

  int n_items = rankings.n_elem;

  double prob = 1;
  while(unranked_items.n_elem > 0) {
    vec available_rankings = rankings(unranked_items);
    int item_to_rank = unranked_items(0);

    vec log_numerator(available_rankings.n_elem);
    for(int ll{}; ll < available_rankings.n_elem; ll++) {
      log_numerator(ll) = -alpha / n_items *
        get_rank_distance(rho(span(item_to_rank)), available_rankings(span(ll)), metric);
    }
    vec sample_probs = normalize_weights(log_numerator);
    if(forward) {

      rankings(span(item_to_rank)) = sample(available_rankings, 1, false, sample_probs);
    }

    int ranking_chosen = as_scalar(find(rankings(item_to_rank) == available_rankings));

    prob *= sample_probs(ranking_chosen);
    if(available_rankings.n_elem <= 1) break;
    unranked_items = unranked_items.subvec(1, available_rankings.n_elem - 1);

    rankings(unranked_items) = setdiff(available_rankings, available_rankings(span(ranking_chosen)));
  }

  return Rcpp::List::create(
    Rcpp::Named("proposal") = rankings,
    Rcpp::Named("probability") = prob
  );
}
