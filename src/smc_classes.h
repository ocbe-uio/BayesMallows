#pragma once

#include <string>
#include "classes.h"
#include "resampler.h"

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

struct LatentParticle {
  LatentParticle(const arma::mat& augmented_data,
                 const arma::uvec& particle_consistent,
                 const unsigned int n_assessors);
  arma::mat augmented_data;
  arma::uvec consistent;
  double aug_acceptance{};
  int aug_count{};
  arma::vec log_proposal_prob;
  double log_inc_wgt{};
};

struct StaticParticle {
  StaticParticle(double alpha, const arma::vec& rho, const unsigned int n_assessors,
                 const std::vector<LatentParticle>& lp);
  ~StaticParticle() = default;

  void prepare_particle_filter(const SMCData& dat);

  std::vector<LatentParticle> particle_filters;
  double alpha;
  arma::vec rho;
  double log_inc_wgt{};
  arma::vec previous_distance;
  double marginal_log_likelihood{};
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

  const unsigned int mcmc_steps;
  const std::string metric;
  const unsigned int n_particles;
  const unsigned int n_particle_filters;
  const int resampling_threshold;
  const std::unique_ptr<RhoProposal> rho_proposal_function;
  const double alpha_prop_sd;
  std::string resampling_method;
};

struct SMCAugmentation {
  SMCAugmentation(
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options
    );
  ~SMCAugmentation() = default;

  void run_particle_filter(
      std::vector<StaticParticle>& pvec,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun,
      const std::unique_ptr<Resampler>& resampler,
      size_t time) const;

  void update_missing_ranks(StaticParticle& p, const SMCData& dat,
      const std::unique_ptr<Distance>& distfun) const;

  const int max_topological_sorts;

private:
  const std::unique_ptr<PartialProposal> partial_aug_prop;
  const std::unique_ptr<PairwiseProposal> pairwise_aug_prop;
  const unsigned int latent_sampling_lag;
  StaticParticle augment_partial(const StaticParticle& pvec, const SMCData& dat) const;
};
