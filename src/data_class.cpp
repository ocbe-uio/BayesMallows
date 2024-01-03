#include "classes.h"
using namespace arma;

triply_nested define_items(const Rcpp::List& data, const std::string type) {
    Rcpp::List input_constraints(data["constraints"]);
  triply_nested items{};
  for(Rcpp::List constraint : input_constraints) {
    doubly_nested user_items{};
    Rcpp::List input = constraint[type];
    for(std::vector<unsigned int> one_input : input) {
      user_items.push_back(one_input);
    }
    items.push_back(user_items);
  }
  return items;
}

Data::Data(
  const Rcpp::List& data
) :
  rankings(Rcpp::as<mat>(data["rankings"]).t()),
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  observation_frequency(data["observation_frequency"]),
  items_above { define_items(data, "items_above") },
  items_below { define_items(data, "items_below") },
  any_missing { !is_finite(rankings) }
  {
    rankings.replace(datum::nan, 0);
  }
