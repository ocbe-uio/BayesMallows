#ifndef PARAMETERUPDATES_H
#define PARAMETERUPDATES_H

#include "RcppArmadillo.h"
#include "classes.h"

struct AlphaRatio{
  AlphaRatio(double proposal, bool accept) :
  proposal {proposal}, accept {accept} {}
  ~AlphaRatio() = default;
  double proposal;
  bool accept;
};

AlphaRatio make_new_alpha(const double& alpha_old, const arma::vec& rho_old,
                          const double& alpha_prop_sd, const std::string& metric,
                          const Rcpp::List& logz_list,
                          const arma::mat& rankings,
                          const arma::vec& observation_frequency,
                          const double& n_items,
                          const Priors& priors);

arma::vec make_new_rho(arma::vec current_rho, const arma::mat& rankings, double alpha_old, int leap_size, std::string metric,
                 arma::vec observation_frequency);

#endif
