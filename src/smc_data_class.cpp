#include "smc_classes.h"
using namespace arma;

SMCData::SMCData(const Rcpp::List& data) :
  Data(data),
  timepoint { Rcpp::as<uvec>(data["timepoint"]) },
  consistent ( data["consistent"]),
  user_ids ( Rcpp::as<Rcpp::Nullable<Rcpp::CharacterVector>>(data["user_ids"]) ) {}

void SMCData::update(const Rcpp::List& new_data) {
  SMCData new_dat{new_data};
  new_rankings = new_dat.rankings;
  num_new_obs = new_dat.n_assessors;
  if(n_assessors == 0) {
    timepoint = new_dat.timepoint;
    rankings = new_dat.rankings;
    missing_indicator = new_dat.missing_indicator;
    observation_frequency = new_dat.observation_frequency;
    consistent = new_dat.consistent;
  } else if (new_dat.n_assessors > 0){
    timepoint = join_cols(timepoint, new_dat.timepoint);
    rankings = join_rows(rankings, new_dat.rankings);
    missing_indicator = join_rows(missing_indicator, new_dat.missing_indicator);
    observation_frequency = join_cols(observation_frequency, new_dat.observation_frequency);
    consistent = join_cols(consistent, new_dat.consistent);
  }
  n_assessors += new_dat.n_assessors;
}
