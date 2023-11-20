#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>
#include "parameterupdates.h"
#include "misc.h"
#include "distances.h"
#include "partitionfuns.h"

struct parameters {
  parameters(Rcpp::List model, Rcpp::List compute_options,
             Rcpp::List priors, Rcpp::List initial_values, int n_items);
  ~parameters() = default;

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

private:
  int kappa_1;
  int kappa_2;
  const int n_items;
  int leap_size;
  int rho_thinning;
  std::string metric;
  double lambda;
  double alpha_prop_sd;
};


#endif
