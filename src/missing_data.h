#pragma once
#include "distances.h"
#include "classes.h"

arma::vec make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    const double& alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<Distance>& pseudo_aug_distance,
    double& log_aug_prob,
    bool pseudo = false);

arma::umat set_up_missing(const Data& dat) noexcept;

arma::mat initialize_missing_ranks(
    const arma::mat& rankings,
    const arma::umat& missing_indicator);

struct RankProposal{
  RankProposal() {};
  RankProposal(
    const arma::vec& rankings,
    const double& probability,
    const arma::uvec& mutated_items) :
  rankings { rankings }, probability { probability },
  mutated_items { mutated_items } {}
  ~RankProposal() = default;

  arma::vec rankings{};
  double probability{};
  arma::uvec mutated_items{};
};

RankProposal make_uniform_proposal(
    const arma::vec& ranks,
    const arma::uvec& indicator) noexcept;

RankProposal make_pseudo_proposal(
    arma::vec ranks,
    const arma::uvec& indicator,
    const double& alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun
) noexcept;
