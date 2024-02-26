#pragma once

#include "classes.h"
#include "resampler.h"

struct Particle {
  Particle(double alpha, const arma::vec& rho, const arma::mat& augmented_data,
           const unsigned int n_assessors, const arma::uvec& particle_consistent);
  ~Particle() = default;

  double alpha;
  arma::vec rho;
  arma::mat augmented_data;
  double log_inc_wgt{};
  arma::vec log_aug_prob;
  arma::uvec consistent;
  arma::vec previous_distance;
};

struct SMCData : Data {
  SMCData(const Rcpp::List& data);
  void update(const Rcpp::List& new_data);
  arma::mat new_rankings{};
  unsigned int num_new_obs{};
  arma::uvec timepoint{};
  arma::umat consistent{};
};

struct SMCParameters {
  SMCParameters(
    const Rcpp::List& model_options,
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options);
  ~SMCParameters() = default;

  void update_alpha(
      Particle& p,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun,
      const Priors& priors) const;

  void update_rho(
      Particle& p,
      const SMCData& dat,
      const std::unique_ptr<Distance>& distfun,
      const std::unique_ptr<ProposalDistribution>& prop
  ) const;

  arma::ivec draw_resampling_index(
      const std::vector<Particle>& pvec) const;
  void resample(std::vector<Particle>& pvec) const;

  const unsigned int mcmc_steps;
  const int leap_size;
  const std::string rho_proposal_option;
  const std::string metric;

private:
  const double alpha_prop_sd;
  const std::unique_ptr<Resampler> resampler;
};

struct SMCAugmentation {
  SMCAugmentation(
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options
    );
  ~SMCAugmentation() = default;

  void reweight(
      std::vector<Particle>& pvec,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun) const;

  void augment_partial(std::vector<Particle>& pvec, const SMCData& dat) const;

  void update_missing_ranks(Particle& p, const SMCData& dat,
      const std::unique_ptr<Distance>& distfun) const;

  const std::string aug_method;
  const std::string pseudo_aug_metric;

private:
  const unsigned int latent_sampling_lag;
};
