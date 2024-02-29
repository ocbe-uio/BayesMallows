#pragma once
#include <RcppArmadillo.h>
#include "partition_functions.h"
#include "distances.h"
#include "rank_proposal.h"

using doubly_nested = std::vector<std::vector<unsigned int>>;
using triply_nested = std::vector<doubly_nested>;

struct Data {
  Data(const Rcpp::List& data);
  ~Data() = default;

  arma::mat rankings;
  unsigned int n_assessors;
  const unsigned int n_items;
  arma::vec observation_frequency;
  const triply_nested items_above{};
  const triply_nested items_below{};
  const bool any_missing;
  const bool augpair;
  arma::umat missing_indicator;
};

struct Priors {
  Priors(
    const Rcpp::List& priors
  ) : gamma { priors["gamma"] }, lambda { priors["lambda"] },
  kappa(priors["kappa"]), psi { priors["psi"] } {}
  ~Priors() = default;

  const double gamma;
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

  void update_shape(const Data& dat, const Priors& priors);
  void update_rho(const Data& dat,
                  const arma::uvec& cluster_assignment,
                  const std::unique_ptr<Distance>& distfun,
                  const std::unique_ptr<ProposalDistribution>& prop);

  void update_alpha(
      const Data& dat,
      const std::unique_ptr<Distance>& distfun,
      const std::unique_ptr<PartitionFunction>& pfun,
      const Priors& priors,
      const arma::uvec& current_cluster_assignment);

  arma::mat alpha;
  double alpha_acceptance{};
  arma::vec alpha_old;
  arma::cube rho;
  double rho_acceptance{};
  arma::mat rho_old;
  arma::vec shape_1;
  arma::vec shape_2;
  arma::vec theta;
  double theta_current;
  const unsigned int n_clusters;
  const unsigned int nmc;
  const std::string error_model;
  const int alpha_jump;
  const int leap_size;
  const std::string rho_proposal_option;
  const std::string metric;
  int burnin{};
  size_t t{};
  size_t alpha_index{};
  size_t rho_index{};

private:
  const double alpha_prop_sd;
  const int rho_thinning;
};

struct Clustering {
  Clustering(const Parameters& pars, const Rcpp::List& compute_options,
             const unsigned int n_assessors);
  ~Clustering() = default;

  void update_cluster_probs(const Parameters& pars, const Priors& pris);
  void update_cluster_labels(const int t, const Data& dat,
                             const Parameters& pars,
                             const std::unique_ptr<PartitionFunction>& pfun);
  void update_wcd(const int t);
  void update_dist_mat(const Data& dat, const Parameters& pars,
                       const std::unique_ptr<Distance>& distfun);
  void save_cluster_parameters(size_t t);

  arma::vec current_cluster_probs;
  const bool clustering;
  const unsigned int clus_thinning;
  size_t cluster_assignment_index{};

private:
  const arma::uvec index;
  arma::mat dist_mat;
  const bool include_wcd;
  const bool save_ind_clus;
  int n_cluster_assignments;

public:
  arma::mat cluster_probs;
  arma::umat cluster_assignment;
  arma::uvec current_cluster_assignment;
  arma::mat within_cluster_distance;
};

struct Augmentation {
  Augmentation(Data& dat, const Rcpp::List& compute_options);
  ~Augmentation() = default;

  void augment_pairwise(
      Data& dat,
      const Parameters& pars,
      const Clustering& clus,
      const std::unique_ptr<Distance>& distfun,
      const std::unique_ptr<ProposalDistribution>& prop);

  void update_missing_ranks(
      Data& dat,
      const Clustering& clus,
      const Parameters& pars,
      const std::unique_ptr<Distance>& distfun
  );
  void save_augmented_data(const Data& dat, const Parameters& pars);
  const bool save_aug;
  const unsigned int aug_thinning;
  arma::cube augmented_data;
  const unsigned int swap_leap;
  size_t aug_index{};

private:
  const std::string aug_method;
  const std::string pseudo_aug_metric;
  const std::unique_ptr<Distance> pseudo_aug_distance;
  arma::vec log_aug_prob;
};

