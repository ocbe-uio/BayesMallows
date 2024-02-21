#pragma once
#include "distances.h"
#include "classes.h"

arma::vec make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    double alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<Distance>& pseudo_aug_distance,
    double& log_aug_prob);

arma::umat set_up_missing(const Data& dat) noexcept;

arma::mat initialize_missing_ranks(
    arma::mat rankings,
    const arma::umat& missing_indicator);

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
};

RankProposal make_uniform_proposal(
    const arma::vec& ranks,
    const arma::uvec& indicator) noexcept;

RankProposal make_pseudo_proposal(
    arma::vec ranks,
    const arma::uvec& indicator,
    double alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun
) noexcept;
