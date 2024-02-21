#include <memory>
#include <utility>
#include "rank_proposal.h"
#include "distances.h"
#include "missing_data.h"

using namespace arma;

std::unique_ptr<ProposalDistribution> choose_rank_proposal(
    const std::string& rho_proposal, int leap_size,
    const std::unique_ptr<Distance>& distfun) {
  if(rho_proposal == "ls") {
    return std::make_unique<LeapAndShift>(leap_size, distfun);
  } else if(rho_proposal == "swap") {
    return std::make_unique<Swap>(leap_size, distfun);
  } else {
    Rcpp::stop("Unknown proposal distribution.");
  }
}

ProposalDistribution::ProposalDistribution(
  int leap_size, const std::unique_ptr<Distance>& distfun) :
  leap_size { leap_size }, distfun(distfun) {}

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
    const vec& current_rank) {
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

std::pair<unsigned int, unsigned int> Swap::sample(
    const vec& current_rank) {
  int n_items = current_rank.n_elem;
  ivec l = Rcpp::sample(leap_size, 1);
  ivec a = Rcpp::sample(n_items - l(0), 1);
  int u = a(0);

  unsigned int ind1 = as_scalar(find(current_rank == u));
  unsigned int ind2 = as_scalar(find(current_rank == (u + l(0))));
  return std::make_pair(ind1, ind2);
}

RankProposal Swap::propose(const vec& current_rank) {
  auto inds = sample(current_rank);
  RankProposal rp{current_rank};
  std::swap(rp.rankings(inds.first), rp.rankings(inds.second));
  rp.mutated_items = {inds.first, inds.second};
  distfun->update_leap_and_shift_indices(rp.mutated_items, current_rank.n_elem);

  return rp;
}

RankProposal Swap::propose(
    const arma::vec& current_rank,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below) {
  auto inds = sample(current_rank);
  RankProposal ret{current_rank};
  std::swap(ret.rankings(inds.first), ret.rankings(inds.second));

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

    ret.g_diff += count_error_diff(inds.first) + count_error_diff(inds.second);
    return ret;
}
