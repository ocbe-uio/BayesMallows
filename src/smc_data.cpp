#include "classes.h"
using namespace arma;

SMCData::SMCData(
  const Rcpp::List& data,
  const Rcpp::List& new_data
) : Data(data),
new_rankings { Rcpp::as<mat>(new_data["rankings"]).t() },
num_new_obs { new_rankings.n_cols },
consistent (data["consistent"]) {}
