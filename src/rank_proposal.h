#pragma once
#include <RcppArmadillo.h>

struct Distance;
struct RankProposal;

struct ProposalDistribution {
  ProposalDistribution(int leap_size);
  virtual ~ProposalDistribution() = default;
  virtual RankProposal propose(
      const arma::vec& current_rank,
      const std::unique_ptr<Distance>& distfun) = 0;

  int leap_size;
};

std::unique_ptr<ProposalDistribution> choose_rank_proposal(
    const std::string& rho_proposal, int leap_size);

struct LeapAndShift : ProposalDistribution {
  using ProposalDistribution::ProposalDistribution;

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
    const arma::vec& current_rank,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below);

};

struct Swap : ProposalDistribution {
  using ProposalDistribution::ProposalDistribution;
};
