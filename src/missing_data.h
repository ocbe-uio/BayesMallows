#pragma once
#include "distances.h"
#include "rank_proposal.h"

arma::vec make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    double alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartialProposal>& propfun,
    double& log_aug_prob);

arma::mat initialize_missing_ranks(
    arma::mat rankings,
    const arma::umat& missing_indicator);
