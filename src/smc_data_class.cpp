#include "smc_classes.h"
using namespace arma;

SMCData::SMCData(
  const Rcpp::List& data,
  const Rcpp::List& new_data
) : Data(data),
new_rankings { Rcpp::as<mat>(new_data["rankings"]).t() },
num_new_obs { new_rankings.n_cols } {}

void SMCData::update_data(const Particle& p) {
  if(!any_missing) return;
  rankings = p.augmented_data;
}
