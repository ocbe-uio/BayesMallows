#pragma once

#include <string>
#include "classes.h"
#include "resampler.h"

struct StaticParticle {
  StaticParticle(double alpha, const arma::vec& rho, const arma::mat& augmented_data,
           const unsigned int n_assessors, const arma::uvec& particle_consistent);
  ~StaticParticle() = default;

  double alpha;
  arma::vec rho;
  arma::mat augmented_data;
  double log_inc_wgt{};
  arma::vec log_aug_prob;
  arma::uvec consistent;
  arma::vec previous_distance;
  double alpha_acceptance{};
  double rho_acceptance{};
  double aug_acceptance{};
  int aug_count{};
};

struct SMCData : Data {
  SMCData(const Rcpp::List& data);
  void update(const Rcpp::List& new_data);
  Rcpp::List wrapup();
  arma::mat new_rankings{};
  unsigned int num_new_obs{};
  arma::uvec timepoint{};
  arma::umat consistent{};
  Rcpp::IntegerVector user_ids{};
  Rcpp::IntegerVector updated_match{};
  arma::imat preferences;
};

struct SMCParameters {
  SMCParameters(
    const Rcpp::List& model_options,
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options);
  ~SMCParameters() = default;

  void update_alpha(
      StaticParticle& p,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun,
      const Priors& priors) const;

  void update_rho(
      StaticParticle& p,
      const SMCData& dat,
      const std::unique_ptr<Distance>& distfun
  ) const;

  void resample(std::vector<StaticParticle>& pvec);
  const unsigned int mcmc_steps;
  const std::string metric;
  const unsigned int n_particles;
  bool rejuvenate{false};

private:
  const std::unique_ptr<RhoProposal> rho_proposal_function;
  const double alpha_prop_sd;
  const std::unique_ptr<Resampler> resampler;
  const int resampling_threshold;
};

struct SMCAugmentation {
  SMCAugmentation(
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options
    );
  ~SMCAugmentation() = default;

  void reweight(
      std::vector<StaticParticle>& pvec,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun) const;

  void update_missing_ranks(StaticParticle& p, const SMCData& dat,
      const std::unique_ptr<Distance>& distfun) const;

  const int max_topological_sorts;

private:
  const std::unique_ptr<PartialProposal> partial_aug_prop;
  const std::unique_ptr<PairwiseProposal> pairwise_aug_prop;
  const unsigned int latent_sampling_lag;
  std::vector<StaticParticle> augment_partial(
      const std::vector<StaticParticle>& pvec, const SMCData& dat) const;
};
