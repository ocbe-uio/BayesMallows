#include <vector>
#include <numeric>
#include "smc_classes.h"
#include "missing_data.h"

using namespace arma;

Particle::Particle(
  double alpha, const vec& rho, const mat& augmented_data,
  const unsigned int n_assessors, const uvec& particle_consistent) :
  alpha (alpha), rho (rho), augmented_data (augmented_data),
  log_aug_prob (zeros(n_assessors)),
  consistent(particle_consistent),
  previous_distance(zeros(n_assessors)){}

std::vector<Particle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    const SMCData& dat
) {
  vec alpha_samples(initial_values["alpha_init"]);
  mat rho_samples(initial_values["rho_init"]);
  Rcpp::Nullable<cube> aug_init(initial_values["aug_init"]);
  if(rho_samples.n_rows != dat.n_items) {
    Rcpp::stop("Wrong format for initial values for rho.");
  }

  std::vector<Particle> pvec;
  pvec.reserve(n_particles);

  for(size_t i{}; i < n_particles; i++) {
    uvec particle_consistent;
    mat augmented_data;
    if(dat.any_missing || dat.augpair) {
      if(aug_init.isNotNull()) {
        particle_consistent = uvec(dat.n_assessors - dat.num_new_obs, fill::ones);
        augmented_data = Rcpp::as<cube>(aug_init).slice(i);
      } else {
        augmented_data = initialize_missing_ranks(dat.rankings, dat.missing_indicator);
      }
    }

    pvec.emplace_back(
      Particle(alpha_samples(i), rho_samples.col(i), augmented_data,
               dat.n_assessors, particle_consistent));
  }

  return pvec;
}

std::vector<Particle> augment_particles(
    const std::vector<Particle>& pvec_init,
    const SMCData& dat
) {
  auto pvec = pvec_init;
  for(size_t i{}; i < pvec.size(); i++) {
    pvec[i].alpha_acceptance = 0;
    pvec[i].rho_acceptance = 0;
    uvec particle_consistent;
    if(dat.any_missing) {
      particle_consistent = uvec(dat.n_assessors - dat.num_new_obs, fill::ones);

      for(auto index : dat.updated_match) {
        vec to_compare = dat.rankings.col(index);
        uvec comparison_inds = find(to_compare > 0);
        vec augmented = pvec[i].augmented_data(span::all, span(index));

        particle_consistent(index) =
          all(to_compare(comparison_inds) == augmented(comparison_inds));
      }

      if(dat.num_new_obs > 0) {
        mat tmp = initialize_missing_ranks(
          dat.new_rankings,
          dat.missing_indicator(
            span::all,
            span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)));

        pvec[i].augmented_data.resize(dat.n_items, dat.rankings.n_cols);
        pvec[i].augmented_data(
            span::all,
            span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)
        ) = tmp;
        pvec[i].log_aug_prob.resize(dat.rankings.n_cols);
      }
    } else if (dat.augpair) {
      pvec[i].augmented_data.resize(dat.n_items, dat.rankings.n_cols);
      pvec[i].augmented_data(
          span::all,
          span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)
      ) = dat.rankings(span::all, span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1));
      pvec[i].log_aug_prob.resize(dat.rankings.n_cols);
    }
  }

  return pvec;
}

cube wrapup_rho(const std::vector<std::vector<Particle>>& particle_vectors) {
  cube rho_samples(particle_vectors[0][0].rho.size(),
                   particle_vectors[0].size(),
                   particle_vectors.size());
  for(size_t i{}; i < particle_vectors.size(); i++) {
    for(size_t j{}; j < particle_vectors[i].size(); j++) {
      rho_samples(span::all, span(j), span(i)) = particle_vectors[i][j].rho;
    }
  }

  return rho_samples;
}

mat wrapup_alpha(const std::vector<std::vector<Particle>>& particle_vectors) {
  mat alpha_samples(particle_vectors[0].size(), particle_vectors.size());
  for(size_t i{}; i < particle_vectors.size(); i++) {
    for(size_t j{}; j < particle_vectors[i].size(); j++) {
      alpha_samples(j, i) = particle_vectors[i][j].alpha;
    }
  }

  return alpha_samples;
}

cube wrapup_augmented_data(const std::vector<Particle>& pvec) {
  cube augmented_data;
  if(!pvec[0].augmented_data.is_empty()) {
    augmented_data.set_size(pvec[0].augmented_data.n_rows,
                            pvec[0].augmented_data.n_cols,
                            pvec.size());
    for(size_t i{}; i < pvec.size(); i++) {
      augmented_data.slice(i) = pvec[i].augmented_data;
    }
  }
  return augmented_data;
}

Rcpp::List compute_particle_acceptance(
    const std::vector<std::vector<Particle>>& particle_vectors, int mcmc_steps) {
  vec alpha_acceptance(particle_vectors.size());
  vec rho_acceptance(particle_vectors.size());
  vec aug_acceptance(particle_vectors.size());
  for(size_t i{}; i < particle_vectors.size(); i++) {
    alpha_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const Particle& p) {
          return accumulator + p.alpha_acceptance;
          });
    alpha_acceptance[i] /= particle_vectors[i].size() * mcmc_steps;
    rho_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const Particle& p) {
          return accumulator + p.rho_acceptance;
        });
    rho_acceptance[i] /= particle_vectors[i].size() * mcmc_steps;
    aug_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const Particle& p) {
          return accumulator + p.aug_acceptance / p.aug_count;
        });
    aug_acceptance[i] /= particle_vectors[i].size();
  }
  return Rcpp::List::create(
    Rcpp::Named("alpha_acceptance") = alpha_acceptance,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("aug_acceptance") = aug_acceptance
  );
}
