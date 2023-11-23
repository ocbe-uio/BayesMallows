#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "misc.h"
#include "distances.h"
#include "partitionfuns.h"
#include "missing_data.h"
#include "sample.h"
#include "pairwise_comparisons.h"

template <typename T>
static T verify_positive(const T input) {
  if(input < 0) {
    Rcpp::stop("Positive value required.\n");
  }
  return input;
}


struct Data {
  Data(const Rcpp::List& data, const Rcpp::List& compute_options);
  ~Data() = default;

  arma::mat rankings;
  const Rcpp::List constraints;
  const unsigned int n_assessors;
  const unsigned int n_items;
  const arma::vec observation_frequency;



};


struct Priors {
  Priors(const Rcpp::List& priors);
  ~Priors() = default;

  const double lambda;
  const unsigned int kappa_1;
  const unsigned int kappa_2;
  const unsigned int psi;
  const double theta_error{0.1};
};

struct Parameters {
  Parameters(const Rcpp::List& model,
             const Rcpp::List& compute_options,
             const Rcpp::List& initial_values,
             const unsigned int n_items);
  ~Parameters() = default;

  void update_shape(int t, const Data& dat, const Priors& priors);
  void update_rho(int cluster_index, int t, int& rho_index,
                  const Data& dat);

  void update_alpha(
      int cluster_index,
      int alpha_index,
      const Data& dat,
      const Rcpp::List& logz_list,
      const Priors& priors);

  arma::mat alpha;
  arma::vec alpha_old;
  arma::cube rho;
  arma::mat rho_old;
  arma::vec shape_1;
  arma::vec shape_2;
  arma::vec theta;

  const int n_clusters;
  const int nmc;
  const std::string metric;
  const std::string error_model;

  const int get_alpha_jump() {
    return alpha_jump;
  }

private:

  static std::string verify_metric(const std::string input) {
    bool check = (input.compare("footrule") == 0) ||
    (input.compare("spearman") == 0) ||
    (input.compare("cayley") == 0) ||
    (input.compare("kendall") == 0) ||
    (input.compare("ulam") == 0) ||
    (input.compare("hamming") == 0);
    if(!check) {
      Rcpp::stop("Unknown metric.\n");
    }
    return input;
  }

  static std::string verify_error_model(const std::string input) {
    bool check = (input.compare("none") == 0) ||
      (input.compare("bernoulli") == 0);
    if(!check) {
      Rcpp::stop("Unknown error model.\n");
    }
    return input;
  }

  const int alpha_jump;
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

  const bool clustering;
  const unsigned int clus_thinning;
  const bool include_wcd;
  const bool save_ind_clus;

  void update_cluster_probs(const Parameters& pars, const Priors& pris);
  void update_cluster_labels(const int t, const Data& dat,
                             const Parameters& pars,
                             const Rcpp::List& logz_list);

  void update_wcd(const int t);
  void update_dist_mat(const Data& dat, const Parameters& pars);


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

  void augment_pairwise(Data& dat, const Parameters& pars,
                        const Clustering& clus, const Priors& pris);

  void update_missing_ranks(
      Data& dat,
      const Clustering& clus,
      const Parameters& pars
  );


};


#endif

