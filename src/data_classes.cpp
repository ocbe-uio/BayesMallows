#include "classes.h"
using namespace arma;

Data::Data(
  const Rcpp::List& data
) :
  rankings { Rcpp::as<mat>(data["rankings"]).t() },
  constraints ( Rcpp::as<Rcpp::List>(data["constraints"]) ),
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  observation_frequency { Rcpp::as<vec>(data["observation_frequency"]) } {}

SMCData::SMCData(
  const Rcpp::List& data,
  const Rcpp::List& new_data
) : Data(data),
new_rankings { Rcpp::as<mat>(new_data["rankings"]).t() },
num_new_obs { new_rankings.n_cols },
consistent {Rcpp::as<umat>(new_data["consistent"])} {}
