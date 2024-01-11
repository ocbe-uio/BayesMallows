#pragma once

#include "classes.h"

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
};

struct SMCData : Data {
  SMCData(const Rcpp::List& data, const Rcpp::List& new_data);
  void update_data(const Particle& p);
  arma::mat new_rankings;
  const unsigned int num_new_obs;
};

struct SMCParameters {
  SMCParameters(
    const Rcpp::List& smc_options, const Rcpp::List& compute_options);
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

  Rcpp::IntegerVector draw_resampling_index(const std::vector<Particle>& pvec);
  void resample(std::vector<Particle>& pvec);

  const unsigned int mcmc_steps;

private:
  const double alpha_prop_sd;
  const unsigned int leap_size;
};

struct SMCAugmentation {
  SMCAugmentation(SMCData& dat, const Rcpp::List& compute_options);
  ~SMCAugmentation() = default;

  void reweight(
      std::vector<Particle>& pvec,
      const SMCData& dat,
      const std::unique_ptr<PartitionFunction>& pfun,
      const std::unique_ptr<Distance>& distfun);

  void augment_partial(
      std::vector<Particle>& pvec,
      const SMCData& dat);

  void update_missing_ranks(Particle& p, const SMCData& dat,
      const std::unique_ptr<Distance>& distfun);

  void resample(const arma::uvec& index, const SMCData& dat);
  const arma::umat missing_indicator;

  const std::string aug_method;
  const std::string pseudo_aug_metric;
};
