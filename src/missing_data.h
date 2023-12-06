#pragma once

arma::vec propose_augmentation(
    const arma::vec& ranks,
    const arma::uvec& indicator);

arma::vec make_new_augmentation(
    const arma::vec& rankings,
    const arma::uvec& missing_indicator,
    const double& alpha,
    const arma::vec& rho,
    const std::string& metric,
    bool pseudo = false);

void set_up_missing(
    arma::mat& rankings,
    arma::umat& missing_indicator);

void initialize_missing_ranks(
    arma::mat& rankings,
    const arma::umat& missing_indicator);

struct PseudoProposal{
  PseudoProposal(
    const arma::vec& rankings,
    const double& probability) :
  rankings { rankings }, probability { probability } {}
  ~PseudoProposal() = default;

  arma::vec rankings;
  double probability;
};

PseudoProposal make_pseudo_proposal(
    arma::uvec unranked_items,
    arma::vec rankings,
    const double& alpha,
    const arma::vec& rho,
    const std::string metric,
    const bool forward = true
);
