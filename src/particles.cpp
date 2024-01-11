#include "smc_classes.h"
#include "missing_data.h"

using namespace arma;

Particle::Particle(
  double alpha, const vec& rho, const mat& augmented_data,
  const unsigned int n_assessors) :
  alpha (alpha), rho (rho), augmented_data (augmented_data),
  log_aug_prob (zeros(n_assessors)) {}

mat initialize_augmented_data(
    const unsigned int particle_index,
    const SMCData& dat,
    const SMCAugmentation& aug,
    const Rcpp::Nullable<cube>& aug_init
) {
  mat augmented_data;
  if(dat.any_missing){
    augmented_data.set_size(dat.n_items, dat.n_assessors);

    if(aug_init.isNotNull()) {
      augmented_data(
        span::all, span(0, dat.rankings.n_cols - dat.num_new_obs - 1)) =
          Rcpp::as<cube>(aug_init).slice(particle_index);

      if(dat.num_new_obs > 0) {
        augmented_data(
          span::all,
          span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)
        ) = initialize_missing_ranks(
          dat.new_rankings,
          aug.missing_indicator(
            span::all,
            span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1))
        );
      }
    } else {
      augmented_data = initialize_missing_ranks(dat.rankings, aug.missing_indicator);
    }
  }
  return augmented_data;
}

std::vector<Particle> initialize_particles(
    const Rcpp::List& initial_values,
    const Rcpp::List& smc_options,
    const SMCAugmentation& aug,
    const SMCData& dat
) {
  const unsigned int n_particles { smc_options["n_particles"] };
  vec alpha_samples(initial_values["alpha_init"]);
  mat rho_samples(initial_values["rho_init"]);
  Rcpp::Nullable<cube> aug_init(initial_values["aug_init"]);

  std::vector<Particle> pvec;
  pvec.reserve(n_particles);

  for(size_t i{}; i < n_particles; i++) {
    mat augmented_data = initialize_augmented_data(i, dat, aug, aug_init);
    pvec.emplace_back(
      Particle(alpha_samples(i), rho_samples.col(i), augmented_data,
               dat.n_assessors)
    );
  }

  return pvec;
}

mat wrapup_rho(const std::vector<Particle>& pvec, const SMCData& dat) {
  mat rho_samples(dat.n_items, pvec.size());
  for(size_t i{}; i < pvec.size(); i++) rho_samples.col(i) = pvec[i].rho;
  return rho_samples;
}

vec wrapup_alpha(const std::vector<Particle>& pvec) {
  vec alpha_samples(pvec.size());
  for(size_t i{}; i < pvec.size(); i++) alpha_samples(i) = pvec[i].alpha;
  return alpha_samples;
}

cube wrapup_augmented_data(const std::vector<Particle>& pvec, const SMCData& dat) {
  cube augmented_data;
  if(dat.any_missing) {
    augmented_data.set_size(dat.n_items, dat.n_assessors, pvec.size());
    for(size_t i{}; i < pvec.size(); i++) {
      augmented_data.slice(i) = pvec[i].augmented_data;
    }
  }
  return augmented_data;
}
