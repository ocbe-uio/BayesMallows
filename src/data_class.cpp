#include "classes.h"
#include "missing_data.h"
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

arma::umat set_up_missing(const mat& rankings, bool any_missing) noexcept {
  if(!any_missing) return arma::umat{};
  arma::umat missing_indicator = conv_to<umat>::from(rankings);
  missing_indicator.transform( [](int val) { return (val == 0) ? 1 : 0; } );
  return missing_indicator;
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
  any_missing { data["any_missing"] },
  missing_indicator { set_up_missing(rankings, any_missing) }
  {
    if(n_assessors <= 0) Rcpp::stop("Must have at least one observation.");
    rankings.replace(datum::nan, 0);
  }
