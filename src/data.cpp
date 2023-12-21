#include "classes.h"
using namespace arma;

Data::Data(
  const Rcpp::List& data
) :
  rankings(Rcpp::as<mat>(data["rankings"]).t()),
  constraints(data["constraints"]),
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  observation_frequency(data["observation_frequency"]) {
    Rcpp::List input_constraints(data["constraints"]);

    for(Rcpp::List constraint : input_constraints) {
      std::vector<std::vector<unsigned int>> user_items_above{};
      std::vector<std::vector<unsigned int>> user_items_below{};

      Rcpp::List input_above = constraint["items_above"];
      for(std::vector<unsigned int> above : input_above) {
        user_items_above.push_back(above);
      }
      items_above.push_back(user_items_above);

      Rcpp::List input_below = constraint["items_below"];
      for(std::vector<unsigned int> below : input_below) {
        user_items_below.push_back(below);
      }
      items_below.push_back(user_items_below);
    }

  }
