#pragma once
#include <RcppArmadillo.h>
#include "missing_data.h"

struct ProposalDistribution {
  ProposalDistribution() {};
  virtual ~ProposalDistribution() = default;
  virtual RankProposal propose(
      const arma::vec& current_rank,
      const std::unique_ptr<Distance>& distfun) = 0;

};

struct LeapAndShift : ProposalDistribution {
  LeapAndShift(int leap_size);

  int find_lower_limit(int item, const arma::uvec& items_above_item,
                       const arma::vec& current_ranking);

  int find_upper_limit(int item, const arma::uvec& items_below_item,
                       const arma::vec& current_ranking);

  RankProposal shift(
      const RankProposal& rp_in, const arma::vec& current_rank, int u);

  RankProposal propose(
      const arma::vec& current_rank,
      const std::unique_ptr<Distance>& distfun);
  RankProposal propose(
    const arma::vec& current_rank, const doubly_nested& items_above,
    const doubly_nested& items_below);

  int leap_size;
};
