#pragma once
#include <RcppArmadillo.h>
#include "distances.h"

using namespace arma;

struct RankProposal {
  RankProposal(const arma::vec& rank_proposal) :
  rank_proposal { rank_proposal } {}
  ~RankProposal() = default;
  arma::vec rank_proposal;
  arma::uvec indices{};
  double prob_backward{};
  double prob_forward{};
};

struct ProposalDistribution {
  ProposalDistribution(int leap_size) : leap_size { leap_size } {};
  virtual ~ProposalDistribution() = default;
  virtual RankProposal propose(
      const arma::vec& current_rank,
      const std::unique_ptr<Distance>& distfun) = 0;

  int leap_size;
};

struct LeapAndShift : ProposalDistribution {

  RankProposal propose(
      const arma::vec& current_rank,
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
    rp.rank_proposal(u) = support(index);

    double support_new = std::min(rp.rank_proposal(u) - 1, leap_size * 1.0) +
      std::min(n_items - rp.rank_proposal(u), leap_size * 1.0);

    double delta_r = rp.rank_proposal(u) - current_rank(u);
    rp.indices = zeros<uvec>(std::abs(delta_r) + 1);
    rp.indices[0] = u;

    if(delta_r > 0){
      for(int k = 1; k <= delta_r; ++k){
        int index = as_scalar(find(current_rank == current_rank(u) + k));
        rp.rank_proposal(index) -= 1;
        rp.indices[k] = index;
      }
    } else if(delta_r < 0) {
      for(int k = (-1); k >= delta_r; --k){
        int index = as_scalar(find(current_rank == current_rank(u) + k));
        rp.rank_proposal(index) += 1;
        rp.indices[-(k)] = index;
      }
    }

    distfun->update_leap_and_shift_indices(rp.indices, n_items);

    if(std::abs(rp.rank_proposal(u) - current_rank(u)) == 1){
      rp.prob_forward = 1.0 / (n_items * support.n_elem) + 1.0 / (n_items * support_new);
      rp.prob_backward = rp.prob_forward;
    } else {
      rp.prob_forward = 1.0 / (n_items * support.n_elem);
      rp.prob_backward = 1.0 / (n_items * support_new);
    }

    return rp;
  }
};
