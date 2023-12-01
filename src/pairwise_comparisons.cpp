#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]



void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const vec& current_ranking) {
  // Find the items which are preferred to the given item
  // Items go from 1, ..., n_items, so must use [item - 1]
  uvec items_above = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[item - 1]);
  uvec items_below = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[item - 1]);

  // If there are any items above, we must find the possible rankings
  if(items_above.n_elem > 0) {
    // Again subtracting 1 because of zero-first indexing
    // Find all the rankings of the items that are preferred to *item*
    vec rankings_above = current_ranking.elem(items_above - 1);
    left_limit = max(rankings_above);
  }

  // If there are any items below, we must find the possible rankings
  if(items_below.n_elem > 0) {
    // Find all the rankings of the items that are disfavored to *item*
    vec rankings_below = current_ranking.elem(items_below - 1);
    right_limit = min(rankings_below);
  }

}

vec propose_pairwise_augmentation(const vec& ranking, const Rcpp::List& assessor_constraints) {
  int n_items = ranking.n_elem;

  // Extract the constraints for this particular assessor
  uvec constrained_items = Rcpp::as<uvec>(assessor_constraints[0]);

  // Sample an integer between 1 and n_items
  int item = randi<int>(distr_param(0, n_items - 1));

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;
  find_pairwise_limits(left_limit, right_limit, item + 1, assessor_constraints, ranking);

  // Now complete the leap step by sampling a new proposal uniformly between
  // left_limit + 1 and right_limit - 1
  int proposed_rank = randi<int>(distr_param(left_limit + 1, right_limit - 1));

  // Assign the proposal to the (item-1)th item
  vec proposal = ranking;
  proposal(item) = proposed_rank;

  uvec indices;

  // Do the shift step
  shift_step(proposal, ranking, item, indices);

  return proposal;
}

vec propose_swap(const vec& ranking, const Rcpp::List& assessor_constraints,
                       int& g_diff, const int& swap_leap) {
  int n_items = ranking.n_elem;

  // Draw a random number, representing an item
  int u = randi<int>(distr_param(1, n_items - swap_leap));

  int ind1 = as_scalar(find(ranking == u));
  int ind2 = as_scalar(find(ranking == (u + swap_leap)));

  vec proposal = ranking;
  proposal(ind1) = ranking(ind2);
  proposal(ind2) = ranking(ind1);

  // First consider the first item that was switched
  uvec items_above = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[ind1]);
  uvec items_below = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[ind1]);

  for(size_t j = 0; j < items_above.n_elem; ++j){
    g_diff += (proposal(items_above[j] - 1) > proposal(ind1)) - (ranking(items_above[j] - 1) > ranking(ind1));
  }
  for(size_t j = 0; j < items_below.n_elem; ++j){
    g_diff += (proposal(items_below[j] - 1) < proposal(ind1)) - (ranking(items_below[j] - 1) < ranking(ind1));
  }

  // Now consider the second item
  items_above = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[ind2]);
  items_below = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[ind2]);

  for(size_t j = 0; j < items_above.n_elem; ++j){
    g_diff += (proposal(items_above[j] - 1) > proposal(ind1)) - (ranking(items_above[j] - 1) > ranking(ind1));
  }
  for(size_t j = 0; j < items_below.n_elem; ++j){
    g_diff += (proposal(items_below[j] - 1) < proposal(ind1)) - (ranking(items_below[j] - 1) < ranking(ind1));
  }
  return proposal;
}


