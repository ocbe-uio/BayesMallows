#include <RcppArmadillo.h>
#include "classes.h"
#include "distances.h"
#include "rank_proposal.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


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


