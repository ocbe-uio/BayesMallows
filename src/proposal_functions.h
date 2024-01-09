#pragma once
#include "classes.h"

struct AlphaRatio{
  AlphaRatio(double proposal, bool accept) :
  proposal { proposal }, accept { accept } {}
  ~AlphaRatio() = default;
  double proposal;
  bool accept;
};

AlphaRatio make_new_alpha(
    const double& alpha_old,
    const arma::vec& rho_old,
    const double& alpha_prop_sd,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartitionFunction>& pfun,
    const arma::mat& rankings,
    const arma::vec& observation_frequency,
    const double& n_items,
    const Priors& priors);

arma::vec make_new_rho(
    const arma::vec& current_rho,
    const arma::mat& rankings,
    double alpha_old,
    int leap_size,
    const std::unique_ptr<Distance>& distfun,
    arma::vec observation_frequency);
