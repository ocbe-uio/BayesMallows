#pragma once
#include "distances.h"

arma::vec make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    const double& alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun,
    double& log_aug_prob,
    bool pseudo = false);

void set_up_missing(
    arma::mat& rankings,
    arma::umat& missing_indicator);

void initialize_missing_ranks(
    arma::mat& rankings,
    const arma::umat& missing_indicator);

struct RankProposal{
  RankProposal() {};
  RankProposal(
    const arma::vec& rankings,
    const double& probability) :
  rankings { rankings }, probability { probability } {}
  ~RankProposal() = default;

  arma::vec rankings{};
  double probability{};
};

RankProposal make_uniform_proposal(
    const arma::vec& ranks,
    const arma::uvec& indicator);

RankProposal make_pseudo_proposal(
    arma::uvec unranked_items,
    arma::vec rankings,
    const double& alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun
);
