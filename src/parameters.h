#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "misc.h"
#include "distances.h"
#include "partitionfuns.h"

template <typename T>
static T verify_positive(const T input) {
  if(input < 0) {
    Rcpp::stop("Positive value required.\n");
  }
  return input;
}

struct Data {
  Data(const Rcpp::List& data);
  ~Data() = default;

  arma::mat rankings;
  const unsigned int n_assessors;
  const unsigned int n_items;

};

struct Priors {
  Priors(const Rcpp::List& priors);
  ~Priors() = default;

  const double lambda;
  const unsigned int kappa_1;
  const unsigned int kappa_2;
};

struct Parameters {
  Parameters(const Rcpp::List& model,
             const Rcpp::List& compute_options,
             const Rcpp::List& initial_values,
             const unsigned int n_items);
  ~Parameters() = default;

  void update_shape(int t, const arma::mat& rankings,
                    const Rcpp::List& constraints,
                    const Priors& priors);
  void update_rho(int cluster_index, int t, int& rho_index,
                  const arma::mat& rankings,
                  const arma::vec& observation_frequency);

  void update_alpha(
      int cluster_index,
      int alpha_index,
      const arma::mat& rankings,
      const arma::vec& observation_frequency,
      const Rcpp::List& logz_list,
      const Priors& priors);

  arma::mat alpha;
  arma::vec alpha_old;
  arma::cube rho;
  arma::mat rho_old;
  arma::vec shape_1;
  arma::vec shape_2;
  arma::vec theta;

  const int get_alpha_jump() {
    return alpha_jump;
  }
  const int get_nmc() {
    return nmc;
  }
  const int get_n_clusters() {
    return n_clusters;
  }
  const std::string get_error_model() {
    return error_model;
  }
  const std::string get_metric() {
    return metric;
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
  const std::string error_model;
  const int leap_size;
  const std::string metric;
  const int n_clusters;
  const int nmc;
  const int rho_thinning;

};


#endif
