#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]



void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const uvec& items_above_item,
                          const uvec& items_below_item,
                          const vec& current_ranking) {
  if(items_above_item.size() > 0) {
    vec rankings_above = current_ranking.elem(items_above_item - 1);
    left_limit = max(rankings_above);
  }

  if(items_below_item.size() > 0) {
    vec rankings_below = current_ranking.elem(items_below_item - 1);
    right_limit = min(rankings_below);
  }
}

vec propose_pairwise_augmentation(
    const vec& ranking, const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below) {
  int n_items = ranking.n_elem;

  Rcpp::IntegerVector a = Rcpp::sample(n_items, 1) - 1;
  int item = a(0);

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;
  find_pairwise_limits(left_limit, right_limit, item, items_above[item],
                       items_below[item], ranking);

  // Now complete the leap step by sampling a new proposal uniformly between
  // left_limit + 1 and right_limit - 1
  Rcpp::IntegerVector b = Rcpp::seq(left_limit + 1, right_limit - 1);
  Rcpp::IntegerVector d = Rcpp::sample(b, 1);
  int proposed_rank = d(0);

  // Assign the proposal to the (item-1)th item
  vec proposal = ranking;
  proposal(item) = proposed_rank;

  uvec indices;

  // Do the shift step
  shift_step(proposal, ranking, item, indices);

  return proposal;
}

int increment_g(
    const vec& ranking, const vec& proposal,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below,
    int ind1, int ind2) {
  uvec items_above_item = items_above[ind2];
  uvec items_below_item = items_below[ind2];
  int result{};
  for(size_t j = 0; j < items_above_item.size(); ++j){
    result += (proposal(items_above_item(j) - 1) > proposal(ind1)) -
      (ranking(items_above_item(j) - 1) > ranking(ind1));
  }
  for(size_t j = 0; j < items_below_item.size(); ++j){
    result += (proposal(items_below_item(j) - 1) < proposal(ind1)) -
      (ranking(items_below_item(j) - 1) < ranking(ind1));
  }
  return result;
}

vec propose_swap(
    const vec& ranking,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below,
    int& g_diff, const int& swap_leap) {
  int n_items = ranking.n_elem;

  // Draw a random number, representing an item
  Rcpp::IntegerVector a = Rcpp::sample(n_items - swap_leap, 1);
  int u = a(0);

  int ind1 = as_scalar(find(ranking == u));
  int ind2 = as_scalar(find(ranking == (u + swap_leap)));
  vec proposal = ranking;
  proposal(ind1) = ranking(ind2);
  proposal(ind2) = ranking(ind1);

  g_diff += increment_g(ranking, proposal, items_above,
                        items_below, ind1, ind1);
  g_diff += increment_g(ranking, proposal, items_above,
                        items_below, ind1, ind2);

  return proposal;
}


