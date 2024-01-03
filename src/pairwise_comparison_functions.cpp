#include <RcppArmadillo.h>
#include "classes.h"
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
    const vec& ranking, const doubly_nested& items_above,
    const doubly_nested& items_below) {
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

vec propose_swap(
    const vec& ranking,
    const doubly_nested& items_above,
    const doubly_nested& items_below,
    int& g_diff, const int& swap_leap) {
  int n_items = ranking.n_elem;
  Rcpp::IntegerVector l = Rcpp::sample(swap_leap, 1);
  Rcpp::IntegerVector a = Rcpp::sample(n_items - l(0), 1);
  int u = a(0);

  int ind1 = as_scalar(find(ranking == u));
  int ind2 = as_scalar(find(ranking == (u + l(0))));
  vec proposal = ranking;
  proposal(ind1) = ranking(ind2);
  proposal(ind2) = ranking(ind1);

  auto count_error_diff =
    [&items_above, &items_below, &ranking, &proposal](int ind) {
    int result{};
    for(const auto item_above_first : items_above[ind]) {
      result += (proposal(item_above_first - 1) > proposal(ind)) -
        (ranking(item_above_first - 1) > ranking(ind));
    }
    for(const auto item_below_first : items_below[ind]) {
      result += (proposal(item_below_first - 1) < proposal(ind)) -
        (ranking(item_below_first - 1) < ranking(ind));
    }
    return result;
  };

  g_diff += count_error_diff(ind1) + count_error_diff(ind2);
  return proposal;
}


