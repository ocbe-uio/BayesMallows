#pragma once
#include <RcppArmadillo.h>

struct Distance;
struct RankProposal;

struct RhoProposal {
  RhoProposal(int leap_size);
  virtual ~RhoProposal() = default;
  virtual RankProposal propose(const arma::vec& current_rank) = 0;
  const int leap_size;
};

struct RhoLeapAndShift : RhoProposal {
  using RhoProposal::RhoProposal;
  RankProposal propose(const arma::vec& current_rank) override;
};

struct RhoSwap : RhoProposal {
  using RhoProposal::RhoProposal;
  RankProposal propose(const arma::vec& current_rank) override;
};


struct PairwiseProposal {
  PairwiseProposal();
  virtual ~PairwiseProposal() = default;
  virtual RankProposal propose(
      const arma::vec& current_rank,
      const std::vector<std::vector<unsigned int>>& items_above,
      const std::vector<std::vector<unsigned int>>& items_below) = 0;
};

struct PairwiseLeapAndShift : PairwiseProposal {
  PairwiseLeapAndShift();
  RankProposal propose(
      const arma::vec& current_rank,
      const std::vector<std::vector<unsigned int>>& items_above,
      const std::vector<std::vector<unsigned int>>& items_below) override;
};

struct PairwiseSwap : PairwiseProposal {
  PairwiseSwap(int leap_size);
  RankProposal propose(
      const arma::vec& current_rank,
      const std::vector<std::vector<unsigned int>>& items_above,
      const std::vector<std::vector<unsigned int>>& items_below) override;
  const int leap_size;
};

std::unique_ptr<RhoProposal> choose_rho_proposal(
    const std::string& rho_proposal, int leap_size);

std::unique_ptr<PairwiseProposal> choose_pairwise_proposal(
    const std::string& error_model, unsigned int swap_leap
);
