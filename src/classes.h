#pragma once
#include <RcppArmadillo.h>
#include "partition_functions.h"
#include "distances.h"

struct Data {
  Data(const Rcpp::List& data);
  ~Data() = default;

  arma::mat rankings;
  const Rcpp::List constraints;
  const unsigned int n_assessors;
  const unsigned int n_items;
  const arma::vec observation_frequency;
};

struct Priors {
  Priors(
    const Rcpp::List& priors
  ) : lambda { priors["lambda"] }, kappa(priors["kappa"]),
  psi { priors["psi"] } {}
  ~Priors() = default;

  const double lambda;
  const arma::ivec kappa;
  const unsigned int psi;
};

struct Parameters {
  Parameters(
    const Rcpp::List& model_options,
    const Rcpp::List& compute_options,
    const Rcpp::List& initial_values,
    const unsigned int n_items);
  ~Parameters() = default;

  void update_shape(int t, const Data& dat, const Priors& priors);
  void update_rho(int t, int& rho_index, const Data& dat,
                  const arma::uvec& cluster_assignment,
                  const std::unique_ptr<Distance>& distfun);

  void update_alpha(
      int alpha_index,
      const Data& dat,
      const std::unique_ptr<Distance>& distfun,
      const std::unique_ptr<PartitionFunction>& pfun,
      const Priors& priors,
      const arma::uvec& current_cluster_assignment);

  arma::mat alpha;
  arma::vec alpha_old;
  arma::cube rho;
  arma::mat rho_old;
  arma::vec shape_1;
  arma::vec shape_2;
  arma::vec theta;
  const unsigned int n_clusters;
  const unsigned int nmc;
  const std::string error_model;
  const int alpha_jump;
  const arma::uvec element_indices;

private:
  const double alpha_prop_sd;
  const int leap_size;
  const int rho_thinning;
};

struct Clustering {
  Clustering(const Parameters& pars, const Rcpp::List& compute_options,
             const unsigned int n_assessors);
  ~Clustering() = default;

  arma::mat cluster_probs;
  arma::vec current_cluster_probs;
  arma::umat cluster_assignment;
  arma::uvec current_cluster_assignment;
  arma::mat within_cluster_distance;
  arma::mat dist_mat;
  const arma::uvec index;

  const bool clustering;
  const unsigned int clus_thinning;
  const bool include_wcd;
  const bool save_ind_clus;

  void update_cluster_probs(const Parameters& pars, const Priors& pris);
  void update_cluster_labels(const int t, const Data& dat,
                             const Parameters& pars,
                             const std::unique_ptr<PartitionFunction>& pfun);

  void update_wcd(const int t);
  void update_dist_mat(const Data& dat, const Parameters& pars,
                       const std::unique_ptr<Distance>& distfun);
};

struct Augmentation {
  Augmentation(Data& dat, const Rcpp::List& compute_options);
  ~Augmentation() = default;

  const bool augpair;
  const bool any_missing;
  const bool save_aug;
  const unsigned int aug_thinning;
  const unsigned int swap_leap;

  arma::umat missing_indicator{};
  arma::cube augmented_data{};

  void augment_pairwise(
      const unsigned int t,
      Data& dat,
      const Parameters& pars,
      const Clustering& clus,
      const std::unique_ptr<Distance>& distfun);

  void update_missing_ranks(
      Data& dat,
      const Clustering& clus,
      const Parameters& pars,
      const std::unique_ptr<Distance>& distfun
  );
};

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
  void resample(const arma::uvec& index);
};

struct SMCAugmentation {
  SMCAugmentation(
    SMCData& dat,
    const Rcpp::List& smc_options,
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

  void resample(const arma::uvec& index);
  const std::string aug_method;
  arma::mat log_aug_prob;
  bool any_missing;
  arma::cube augmented_data;
  arma::umat missing_indicator;
  Rcpp::Nullable<arma::cube> aug_init;
};
