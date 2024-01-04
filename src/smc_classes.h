#pragma once

#include "classes.h"

struct SMCData : Data {
  SMCData(const Rcpp::List& data, const Rcpp::List& new_data);
  arma::mat new_rankings;
  const unsigned int num_new_obs;
  const arma::umat consistent;
};

struct SMCParameters {
  SMCParameters(
    const Rcpp::List& model_options,
    const Rcpp::List& smc_options,
    const Rcpp::List& compute_options,
    const Rcpp::List& initial_values);
  ~SMCParameters() = default;

  void update_alpha(
      const unsigned int particle_index,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun,
      const Priors& priors);

  void update_rho(
      const unsigned int particle_index,
      const SMCData& dat,
      const std::unique_ptr<Distance>& distfun
  );

  const unsigned int n_particles;
  const unsigned int mcmc_steps;
  arma::vec alpha_samples;
  arma::mat rho_samples;
  const double alpha_prop_sd;
  const unsigned int leap_size;
  arma::vec log_inc_wgt;
  arma::uvec draw_resampling_index();
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
      const SMCData& dat,
      const std::unique_ptr<Distance>& distfun);

  void update_data(
      const unsigned int particle_index,
      SMCData& dat);

  void update_missing_ranks(
      const unsigned int particle_index,
      const SMCData& dat,
      const SMCParameters& pars,
      const std::unique_ptr<Distance>& distfun);

  void resample(const arma::uvec& index, const SMCData& dat);
  const std::string aug_method;
  const std::string pseudo_aug_metric;
  const std::unique_ptr<Distance> pseudo_aug_distance;
  const Rcpp::Nullable<arma::cube> aug_init;
  const arma::umat missing_indicator;
  arma::mat log_aug_prob;
  arma::cube augmented_data;
};
