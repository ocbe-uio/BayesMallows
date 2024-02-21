#include <RcppArmadillo.h>
#include "classes.h"
#include "distances.h"
#include "rank_proposal.h"
#include "missing_data.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


RankProposal propose_swap(
    const vec& current_rank,
    const doubly_nested& items_above,
    const doubly_nested& items_below,
    int swap_leap) {
  int n_items = current_rank.n_elem;
  ivec l = Rcpp::sample(swap_leap, 1);
  ivec a = Rcpp::sample(n_items - l(0), 1);
  int u = a(0);

  int ind1 = as_scalar(find(current_rank == u));
  int ind2 = as_scalar(find(current_rank == (u + l(0))));
  RankProposal ret{current_rank};
  std::swap(ret.rankings(ind1), ret.rankings(ind2));

  auto count_error_diff =
    [&items_above, &items_below, &current_rank, &ret](int ind) {
    int result{};
    for(const auto item_above_first : items_above[ind]) {
      result += (ret.rankings(item_above_first - 1) > ret.rankings(ind)) -
        (current_rank(item_above_first - 1) > current_rank(ind));
    }
    for(const auto item_below_first : items_below[ind]) {
      result += (ret.rankings(item_below_first - 1) < ret.rankings(ind)) -
        (current_rank(item_below_first - 1) < current_rank(ind));
    }
    return result;
  };

  ret.g_diff += count_error_diff(ind1) + count_error_diff(ind2);
  return ret;
}


