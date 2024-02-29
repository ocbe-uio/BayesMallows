#pragma once
#include <RcppArmadillo.h>

struct Distance;
struct RankProposal;

struct ProposalDistribution {
  ProposalDistribution(int leap_size);
  virtual ~ProposalDistribution() = default;
  virtual RankProposal propose(const arma::vec& current_rank) = 0;
  virtual RankProposal propose(
    const arma::vec& current_rank,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below) = 0;

  const int leap_size;
};

std::unique_ptr<ProposalDistribution> choose_rank_proposal(
    const std::string& rho_proposal, int leap_size);

std::unique_ptr<ProposalDistribution> choose_pairwise_proposal(
  const std::string& error_model, unsigned int swap_leap
);

struct LeapAndShift : ProposalDistribution {
  using ProposalDistribution::ProposalDistribution;

  RankProposal shift(
      const RankProposal& rp_in, const arma::vec& current_rank, int u);

  RankProposal propose(const arma::vec& current_rank) override;
  RankProposal propose(
    const arma::vec& current_rank,
    const std::vector<std::vector<unsigned int>>& items_above,
    const std::vector<std::vector<unsigned int>>& items_below) override;

private:
  int find_lower_limit(int item, const arma::uvec& items_above_item,
                       const arma::vec& current_ranking);

  int find_upper_limit(int item, const arma::uvec& items_below_item,
                       const arma::vec& current_ranking);
};

struct Swap : ProposalDistribution {
  using ProposalDistribution::ProposalDistribution;

  RankProposal propose(const arma::vec& current_rank) override;
  RankProposal propose(
      const arma::vec& current_rank,
      const std::vector<std::vector<unsigned int>>& items_above,
      const std::vector<std::vector<unsigned int>>& items_below) override;

private:
  std::pair<unsigned int, unsigned int> sample(const arma::vec& current_rank);
};
