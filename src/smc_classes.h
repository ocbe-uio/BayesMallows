#pragma once

#include "classes.h"

struct SMCData : Data {
  SMCData(const Rcpp::List& data, const Rcpp::List& new_data);
  arma::mat new_rankings;
  const unsigned int num_new_obs;
  const arma::umat consistent;
};

struct Particle {
  Particle(double alpha, const arma::vec& rho);
  ~Particle() = default;

  double alpha;
  arma::vec rho;
  double log_inc_wgt{};
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
      const Priors& priors);

  void update_rho(
      Particle& p,
      const SMCData& dat,
      const std::unique_ptr<Distance>& distfun
  );

  arma::uvec draw_resampling_index(const std::vector<Particle>& pvec);

  const unsigned int mcmc_steps;

private:
  const double alpha_prop_sd;
  const unsigned int leap_size;
};

struct SMCAugmentation {
  SMCAugmentation(
    SMCData& dat,
    const Rcpp::List& compute_options,
    const Rcpp::List& initial_values,
    const unsigned int n_particles);
  ~SMCAugmentation() = default;

  void reweight(
      SMCParameters& pars,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun);

  void augment_partial(
      const SMCParameters& pars,
      const SMCData& dat);

  void update_data(
      const unsigned int particle_index,
      SMCData& dat);

  void update_missing_ranks(
      const unsigned int particle_index,
      const SMCData& dat,
      const SMCParameters& pars,
      const std::unique_ptr<Distance>& distfun);

  void resample(const arma::uvec& index, const SMCData& dat);

  arma::cube augmented_data;

private:
  const arma::umat missing_indicator;
  const std::string aug_method;
  const std::string pseudo_aug_metric;
  const std::unique_ptr<Distance> pseudo_aug_distance;
  const Rcpp::Nullable<arma::cube> aug_init;
  arma::mat log_aug_prob;
};
