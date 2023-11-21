#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "misc.h"
#include "distances.h"
#include "partitionfuns.h"

struct Parameters {
  Parameters(const Rcpp::List& model,
             const Rcpp::List& compute_options,
             const Rcpp::List& priors,
             const Rcpp::List& initial_values,
             const int n_items);
  ~Parameters() = default;

  void update_shape(int t, const arma::mat& rankings,
                    const Rcpp::List& constraints);
  void update_rho(int cluster_index, int t, int& rho_index,
                  const arma::mat& rankings,
                  const arma::vec& observation_frequency);

  void update_alpha(
      int cluster_index,
      int alpha_index,
      const arma::mat& rankings,
      const arma::vec& observation_frequency,
      const Rcpp::List& logz_list);

  arma::mat alpha;
  arma::vec alpha_old;
  arma::cube rho;
  arma::mat rho_old;
  std::string error_model;
  arma::vec theta;
  arma::vec shape_1;
  arma::vec shape_2;

  int get_nmc() {
    return nmc;
  }
  std::string get_metric() {
    return metric;
  }

private:
  static int verify_positive(const int input) {
    if(input < 0) {
      Rcpp::stop("Positive value required.\n");
    }
    return input;
  }
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
  int kappa_1;
  int kappa_2;
  const std::string metric;
  const int n_items;
  const int nmc;
  int leap_size;
  int rho_thinning;
  double lambda;
  double alpha_prop_sd;
};


#endif
