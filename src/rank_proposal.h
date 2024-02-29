#pragma once
#include <RcppArmadillo.h>
#include "distances.h"

struct RankProposal{
  RankProposal() {};
  RankProposal(const arma::vec& rankings) : rankings { rankings } {}
  RankProposal(
    const arma::vec& rankings,
    double prob_forward, double prob_backward,
    const arma::uvec& mutated_items) :
    rankings { rankings }, prob_forward { prob_forward },
    prob_backward { prob_backward }, mutated_items { mutated_items } {}
  ~RankProposal() = default;

  arma::vec rankings{};
  double prob_forward{1};
  double prob_backward{1};
  arma::uvec mutated_items{};
  int g_diff{};
};

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

struct PartialProposal {
  PartialProposal();
  virtual ~PartialProposal() = default;
  virtual RankProposal propose(
    const arma::vec& current_rank, const arma::uvec& indicator,
    double alpha, const arma::vec& rho) = 0;
};

struct PartialUniform : PartialProposal {
  PartialUniform();
  RankProposal propose(
      const arma::vec& current_rank, const arma::uvec& indicator,
      double alpha, const arma::vec& rho) override;
};

struct PartialPseudoProposal : PartialProposal {
  PartialPseudoProposal(const std::string& pseudo_aug_metric);
  RankProposal propose(
      const arma::vec& current_rank, const arma::uvec& indicator,
      double alpha, const arma::vec& rho) override;

  std::unique_ptr<Distance> distfun;
};

std::unique_ptr<RhoProposal> choose_rho_proposal(
    const std::string& rho_proposal, int leap_size);

std::unique_ptr<PairwiseProposal> choose_pairwise_proposal(
    const std::string& error_model, unsigned int swap_leap
);

std::unique_ptr<PartialProposal> choose_partial_proposal(
  const std::string& aug_method, const std::string& pseudo_aug_metric
);
