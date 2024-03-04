#pragma once

#include <string>
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
      const std::unique_ptr<Distance>& distfun
  ) const;

  arma::ivec draw_resampling_index(
      const std::vector<Particle>& pvec) const;
  void resample(std::vector<Particle>& pvec) const;

  const unsigned int mcmc_steps;
  const std::string metric;
  const unsigned int n_particles;

private:
  const std::unique_ptr<RhoProposal> rho_proposal_function;
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
  void augment_pairwise(std::vector<Particle>& pvec, const SMCData& dat) const;

  void update_missing_ranks(Particle& p, const SMCData& dat,
      const std::unique_ptr<Distance>& distfun) const;

private:
  const std::unique_ptr<PartialProposal> partial_aug_prop;
  const std::unique_ptr<PairwiseProposal> pairwise_aug_prop;
  const unsigned int latent_sampling_lag;
};
