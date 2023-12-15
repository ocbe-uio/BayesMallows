#include "classes.h"
using namespace arma;

Data::Data(
  const Rcpp::List& data
) :
  rankings(Rcpp::as<mat>(data["rankings"]).t()),
  constraints(data["constraints"]),
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  observation_frequency(data["observation_frequency"]) {}
