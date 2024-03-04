#include <memory>
#include <utility>
#include "typedefs.h"
#include "setdiff.h"
#include "rank_proposal.h"

using namespace arma;

std::unique_ptr<RhoProposal> choose_rho_proposal(
    const std::string& rho_proposal, int leap_size) {
  if(rho_proposal == "ls") {
    return std::make_unique<RhoLeapAndShift>(leap_size);
  } else if(rho_proposal == "swap") {
    return std::make_unique<RhoSwap>(leap_size);
  } else {
    Rcpp::stop("Unknown proposal distribution.");
  }
}

std::unique_ptr<PairwiseProposal> choose_pairwise_proposal(
    const std::string& error_model, unsigned int swap_leap
) {
  if(error_model == "none"){
    return std::make_unique<PairwiseLeapAndShift>();
  } else if(error_model == "bernoulli"){
    return std::make_unique<PairwiseSwap>(swap_leap);
  } else {
    Rcpp::stop("error_model must be 'none' or 'bernoulli'");
  }
}

std::unique_ptr<PartialProposal> choose_partial_proposal(
    const std::string& aug_method, const std::string& pseudo_aug_metric
) {
  if(aug_method == "uniform") {
    return std::make_unique<PartialUniform>();
  } else if(aug_method == "pseudo") {
    return std::make_unique<PartialPseudoProposal>(pseudo_aug_metric);
  } else {
    Rcpp::stop("augmentation method must be either 'uniform' or 'pseudo'.");
  }
}

int find_upper_limit(int item, const uvec& items_below_item, const vec& current_ranking) {
  if(items_below_item.size() > 0) {
    return min(current_ranking.elem(items_below_item - 1)) - 1;
  } else {
    return current_ranking.size();
  }
}

int find_lower_limit(int item, const uvec& items_above_item,
                     const vec& current_ranking) {
  if(items_above_item.size() > 0) {
    return max(current_ranking.elem(items_above_item - 1)) + 1;
  } else {
    return 1;
  }
}

RankProposal shift(
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

std::pair<unsigned int, unsigned int> sample(
    const vec& current_rank, int leap_size) {
  int n_items = current_rank.n_elem;
  ivec l = Rcpp::sample(leap_size, 1);
  ivec a = Rcpp::sample(n_items - l(0), 1);
  int u = a(0);

  unsigned int ind1 = as_scalar(find(current_rank == u));
  unsigned int ind2 = as_scalar(find(current_rank == (u + l(0))));
  return std::make_pair(ind1, ind2);
}


RhoProposal::RhoProposal(int leap_size) :
  leap_size { leap_size } {}

RankProposal RhoLeapAndShift::propose(
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

  if(std::abs(rp.rankings(u) - current_rank(u)) == 1){
    rp.prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new);
    rp.prob_backward = rp.prob_forward;
  } else {
    rp.prob_forward = 1.0 / (n_items * support.n_elem);
    rp.prob_backward = 1.0 / (n_items * support_new);
  }

  return rp;
}

RankProposal RhoSwap::propose(const vec& current_rank) {
  auto inds = sample(current_rank, leap_size);
  RankProposal rp{current_rank};
  std::swap(rp.rankings(inds.first), rp.rankings(inds.second));
  rp.mutated_items = {inds.first, inds.second};
  return rp;
}

PairwiseProposal::PairwiseProposal() {}
PairwiseLeapAndShift::PairwiseLeapAndShift() {}

RankProposal PairwiseLeapAndShift::propose(
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

PairwiseSwap::PairwiseSwap(int leap_size) : leap_size { leap_size } {}

RankProposal PairwiseSwap::propose(
    const arma::vec& current_rank,
    const doubly_nested& items_above,
    const doubly_nested& items_below) {
  auto inds = sample(current_rank, leap_size);
  RankProposal ret{current_rank};
  std::swap(ret.rankings(inds.first), ret.rankings(inds.second));
  ret.mutated_items = {inds.first, inds.second};

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

PartialProposal::PartialProposal() {}
PartialUniform::PartialUniform() {}

RankProposal PartialUniform::propose(
    const vec& current_rank, const uvec& indicator,
    double alpha, const vec& rho) {
  vec proposal = current_rank;
  uvec missing_inds = find(indicator == 1);
  vec mutable_current_rank = current_rank(missing_inds);
  ivec inds = Rcpp::sample(mutable_current_rank.size(), mutable_current_rank.size()) - 1;
  proposal(missing_inds) = mutable_current_rank(conv_to<uvec>::from(inds));
  return RankProposal(proposal, 1, 1, missing_inds);
}

PartialPseudoProposal::PartialPseudoProposal(
  const std::string& pseudo_aug_metric) :
  distfun { choose_distance_function(pseudo_aug_metric) } {}

std::pair<arma::vec, double> PartialPseudoProposal::propose_pseudo(
    const vec& current_rank, const uvec& unranked_items, const vec& rho,
    double alpha, bool forward) {
  int n_items = current_rank.n_elem;
  vec proposal = current_rank;
  uvec remaining_unranked_items = unranked_items;
  double prob = 1;
  while(remaining_unranked_items.n_elem > 0) {
    vec available_rankings = proposal(remaining_unranked_items);
    int item_to_rank = remaining_unranked_items(0);

    double rho_for_item = rho(item_to_rank);
    vec log_numerator = -alpha / n_items *
      distfun->scalardist(available_rankings, rho_for_item);

    vec sample_probs = normalise(exp(log_numerator), 1);

    if(forward) {
      ivec ans(sample_probs.size());
      R::rmultinom(1, sample_probs.begin(), sample_probs.size(), ans.begin());
      proposal(span(item_to_rank)) = available_rankings(find(ans == 1));
    }

    int ranking_chosen = as_scalar(find(proposal(item_to_rank) == available_rankings));

    prob *= sample_probs(ranking_chosen);
    if(available_rankings.n_elem <= 1) break;
    remaining_unranked_items = remaining_unranked_items.subvec(1, available_rankings.n_elem - 1);

    proposal(remaining_unranked_items) = setdiff(
      available_rankings, available_rankings(span(ranking_chosen)));
  }

  return std::pair<vec, double>(proposal, prob);
}

RankProposal PartialPseudoProposal::propose(
    const vec& current_rank, const uvec& indicator, double alpha, const vec& rho) {
  uvec missing_inds = find(indicator == 1);
  ivec a = Rcpp::sample(missing_inds.size(), missing_inds.size()) - 1;
  uvec unranked_items = missing_inds(conv_to<uvec>::from(a));

  auto forward_proposal = propose_pseudo(
    current_rank, unranked_items, rho, alpha, true);
  auto backward_proposal = propose_pseudo(
    current_rank, unranked_items, rho, alpha, false);

  return RankProposal(
  forward_proposal.first, forward_proposal.second,
  backward_proposal.second, missing_inds);
}
