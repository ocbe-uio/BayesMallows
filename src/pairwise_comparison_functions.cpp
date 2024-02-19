#include <RcppArmadillo.h>
#include "classes.h"
#include "leapandshift.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

int find_lower_limit(int item, const uvec& items_above_item, const vec& current_ranking) {
  if(items_above_item.size() > 0) {
    return max(current_ranking.elem(items_above_item - 1)) + 1;
  } else {
    return 1;
  }
}

int find_upper_limit(int item, const uvec& items_below_item, const vec& current_ranking) {
  if(items_below_item.size() > 0) {
    return min(current_ranking.elem(items_below_item - 1)) - 1;
  } else {
    return current_ranking.size();
  }
}

vec propose_pairwise_augmentation(
    const vec& ranking, const doubly_nested& items_above,
    const doubly_nested& items_below) {
  int n_items = ranking.n_elem;

  ivec a = Rcpp::sample(n_items, 1) - 1;
  int item = a(0);

  int lower_limit = find_lower_limit(item, items_above[item], ranking);
  int upper_limit = find_upper_limit(item, items_below[item], ranking);

  Rcpp::IntegerVector b = Rcpp::seq(lower_limit, upper_limit);
  ivec d = Rcpp::sample(b, 1);
  int proposed_rank = d(0);

  LeapShiftObject ls{ranking};
  ls.rho_proposal(item) = proposed_rank;

  // Do the shift step
  ls.shift_step(ranking, item);
  return ls.rho_proposal;
}

vec propose_swap(
    const vec& ranking,
    const doubly_nested& items_above,
    const doubly_nested& items_below,
    int& g_diff, int swap_leap) {
  int n_items = ranking.n_elem;
  ivec l = Rcpp::sample(swap_leap, 1);
  ivec a = Rcpp::sample(n_items - l(0), 1);
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


