#include <string>
#include "smc_classes.h"
using namespace arma;

SMCData::SMCData(const Rcpp::List& data) :
  Data(data),
  timepoint { Rcpp::as<uvec>(data["timepoint"]) },
  consistent ( data["consistent"]),
  user_ids ( data["user_ids"] ) {}

void SMCData::update(
    const Rcpp::List& new_data) {
  SMCData new_dat{new_data};

  if(new_dat.user_ids.size() > 0) {
    if(augpair) Rcpp::stop("Pairwise prefences not yet supported with user IDs.");
    Rcpp::CharacterVector new_users = Rcpp::setdiff(new_dat.user_ids, user_ids);
    Rcpp::IntegerVector new_match = Rcpp::match(new_users, new_dat.user_ids) - 1;
    Rcpp::IntegerVector new_indices = new_users.size() > 0 ?
    Rcpp::seq(0, new_users.size() - 1) : Rcpp::IntegerVector{};

    Rcpp::CharacterVector updated_users = Rcpp::intersect(new_dat.user_ids, user_ids);
    updated_match = Rcpp::match(updated_users, user_ids) - 1;
    Rcpp::IntegerVector updated_new = Rcpp::match(updated_users, new_dat.user_ids) - 1;
    Rcpp::IntegerVector updated_indices = updated_users.size() > 0 ?
    Rcpp::seq(0, updated_users.size() - 1) : Rcpp::IntegerVector{};

    for(auto index : new_indices) {
      user_ids.push_back(new_users[index]);
      new_rankings = join_rows(new_rankings, new_dat.rankings.col(new_match[index]));
      timepoint = join_cols(timepoint, new_dat.timepoint(span(new_match[index])));
      missing_indicator = join_rows(
        missing_indicator,
        new_dat.any_missing ? new_dat.missing_indicator.col(new_match[index]) :
        uvec(new_dat.n_items, fill::zeros));
      observation_frequency = join_cols(
        observation_frequency, new_dat.observation_frequency(span(new_match[index])));
    }

    for(auto index : updated_indices) {
      rankings.col(updated_match[index]) = new_dat.rankings.col(updated_new[index]);
    }
  } else {
    new_rankings = new_dat.rankings;
    timepoint = join_cols(timepoint, new_dat.timepoint);
    missing_indicator = join_rows(missing_indicator, new_dat.missing_indicator);
    observation_frequency = join_cols(observation_frequency, new_dat.observation_frequency);
    for(auto items_above_for_user : new_dat.items_above) {
      items_above.push_back(items_above_for_user);
    }
    for(auto items_below_for_user : new_dat.items_below) {
      items_below.push_back(items_below_for_user);
    }
  }

  rankings = join_rows(rankings, new_rankings);
  num_new_obs = new_rankings.n_cols;
  n_assessors += num_new_obs;
}

Rcpp::List SMCData::wrapup() {
  return Rcpp::List::create(
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("consistent") = consistent,
    Rcpp::Named("constraints") = Rcpp::List(),
    Rcpp::Named("n_items") = n_items,
    Rcpp::Named("rankings") = rankings.t(),
    Rcpp::Named("user_ids") = user_ids,
    Rcpp::Named("observation_frequency") = observation_frequency,
    Rcpp::Named("timepoint") = timepoint
  );
}
