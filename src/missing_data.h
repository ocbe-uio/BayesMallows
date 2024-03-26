#pragma once
#include "distances.h"
#include "rank_proposal.h"

std::pair<arma::vec, bool> make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    double alpha,
    const arma::vec& rho,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartialProposal>& propfun);

std::pair<arma::vec, bool> make_new_augmentation(
    const arma::vec& rankings,
    double alpha,
    const arma::vec& rho,
    double theta,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PairwiseProposal>& pairwise_aug_prop,
    const doubly_nested& items_above,
    const doubly_nested& items_below, const std::string& error_model);

arma::mat initialize_missing_ranks(
    arma::mat rankings,
    const arma::umat& missing_indicator);

arma::vec initialize_missing_ranks_vec(
    arma::vec rankings, const arma::uvec& missing_indicator);
