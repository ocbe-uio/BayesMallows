#include "rank_proposal.h"
#include "distances.h"

using namespace arma;

LeapAndShift::LeapAndShift(int leap_size) : leap_size { leap_size } {};

int LeapAndShift::find_lower_limit(int item, const uvec& items_above_item,
                     const vec& current_ranking) {
  if(items_above_item.size() > 0) {
    return max(current_ranking.elem(items_above_item - 1)) + 1;
  } else {
    return 1;
  }
}

int LeapAndShift::find_upper_limit(int item, const uvec& items_below_item, const vec& current_ranking) {
  if(items_below_item.size() > 0) {
    return min(current_ranking.elem(items_below_item - 1)) - 1;
  } else {
    return current_ranking.size();
  }
}

RankProposal LeapAndShift::shift(
    const RankProposal& rp_in, const vec& current_rank, int u) {
  RankProposal rp{rp_in};
  double delta_r = rp.rankings(u) - current_rank(u);
  rp.mutated_items = zeros<uvec>(std::abs(delta_r) + 1);
  rp.mutated_items[0] = u;

  if(delta_r > 0){
    for(int k = 1; k <= delta_r; ++k){
      int index = as_scalar(find(current_rank == current_rank(u) + k));
      rp.rankings(index) -= 1;
      rp.mutated_items[k] = index;
    }
  } else if(delta_r < 0) {
    for(int k = (-1); k >= delta_r; --k){
      int index = as_scalar(find(current_rank == current_rank(u) + k));
      rp.rankings(index) += 1;
      rp.mutated_items[-(k)] = index;
    }
  }

  return rp;
}

RankProposal LeapAndShift::propose(
    const vec& current_rank,
    const std::unique_ptr<Distance>& distfun) {
  RankProposal rp{current_rank};
  int n_items = current_rank.n_elem;

  ivec a = Rcpp::sample(n_items, 1) - 1;
  int u = a(0);

  vec support = join_cols(
    regspace(std::max(1.0, current_rank(u) - leap_size), 1, current_rank(u) - 1),
    regspace(current_rank(u) + 1, 1, std::min(n_items * 1.0, current_rank(u) + leap_size)));

  ivec b = Rcpp::sample(support.n_elem, 1) - 1;
  int index = b(0);
  rp.rankings(u) = support(index);

  double support_new = std::min(rp.rankings(u) - 1, leap_size * 1.0) +
    std::min(n_items - rp.rankings(u), leap_size * 1.0);

  rp = shift(rp, current_rank, u);
  distfun->update_leap_and_shift_indices(rp.mutated_items, n_items);

  if(std::abs(rp.rankings(u) - current_rank(u)) == 1){
    rp.prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new);
    rp.prob_backward = rp.prob_forward;
  } else {
    rp.prob_forward = 1.0 / (n_items * support.n_elem);
    rp.prob_backward = 1.0 / (n_items * support_new);
  }

  return rp;
}

RankProposal LeapAndShift::propose(
    const vec& current_rank, const doubly_nested& items_above,
    const doubly_nested& items_below
) {
  int n_items = current_rank.n_elem;

  ivec a = Rcpp::sample(n_items, 1) - 1;
  int item = a(0);

  int lower_limit = find_lower_limit(item, items_above[item], current_rank);
  int upper_limit = find_upper_limit(item, items_below[item], current_rank);

  Rcpp::IntegerVector b = Rcpp::seq(lower_limit, upper_limit);
  ivec d = Rcpp::sample(b, 1);
  int proposed_rank = d(0);

  RankProposal rp{current_rank};
  rp.rankings(item) = proposed_rank;
  rp = shift(rp, current_rank, item);
  return rp;
}
