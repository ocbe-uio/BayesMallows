#include <vector>
#include <numeric>
#include "smc_classes.h"
#include "missing_data.h"
#include "all_topological_sorts.h"

using namespace arma;

LatentParticle::LatentParticle(
  const mat& augmented_data, const uvec& particle_consistent,
  const unsigned int n_assessors) :
  augmented_data (augmented_data),
  consistent(particle_consistent),
  log_aug_prob (zeros(n_assessors)) {}

StaticParticle::StaticParticle(
  double alpha, const vec& rho, const unsigned int n_assessors,
  const std::vector<LatentParticle>& lp) :
  particle_filters (lp),
  alpha (alpha), rho (rho),
  previous_distance(zeros(n_assessors)) {}

std::vector<StaticParticle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    const SMCData& dat
) {
  vec alpha_samples(initial_values["alpha_init"]);
  mat rho_samples(initial_values["rho_init"]);

  if(rho_samples.n_rows != dat.n_items) {
    Rcpp::stop("Wrong format for initial values for rho.");
  }

  std::vector<StaticParticle> pvec;
  pvec.reserve(n_particles);

  Rcpp::Nullable<cube> aug_init(initial_values["aug_init"]);
  for(size_t i{}; i < n_particles; i++) {
    std::vector<LatentParticle> latvec;
    latvec.reserve(1);
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

    latvec.emplace_back(LatentParticle(augmented_data, particle_consistent, dat.n_assessors));
    pvec.emplace_back(
      StaticParticle(
        alpha_samples(i), rho_samples.col(i), dat.n_assessors,
        latvec));
  }

  return pvec;
}

std::vector<StaticParticle> augment_particles(
    const std::vector<StaticParticle>& pvec_init,
    const SMCData& dat, const SMCAugmentation& aug
) {
  auto pvec = pvec_init;

  std::vector<imat> sorts(dat.user_ids.size());
  if(dat.augpair) {
    for(int i{}; i < dat.n_assessors; i++) {
      Rcpp::IntegerVector test = Rcpp::intersect(dat.updated_match, Rcpp::IntegerVector{i});
      if(i >= (dat.n_assessors - dat.num_new_obs) || test.size() > 0) {
        uvec indices = find(dat.preferences.col(0) == dat.user_ids[i]);
        imat prefs = dat.preferences.rows(indices);
        sorts[i] = all_topological_sorts(prefs.cols(1, 2), dat.n_items,
                                         aug.max_topological_sorts);
      }
    }
  }

  for(size_t static_particle_i{}; static_particle_i < pvec.size(); static_particle_i++) {
    pvec[static_particle_i].alpha_acceptance = 0;
    pvec[static_particle_i].rho_acceptance = 0;

    for(size_t latent_particle_j{}; latent_particle_j < pvec[static_particle_i].particle_filters.size(); latent_particle_j++) {
      if(dat.any_missing || dat.augpair) {
        pvec[static_particle_i].particle_filters[latent_particle_j].consistent = ones<uvec>(dat.n_assessors - dat.num_new_obs);
      }

      if(dat.any_missing) {
        for(auto index : dat.updated_match) {
          vec to_compare = dat.rankings.col(index);
          uvec comparison_inds = find(to_compare > 0);
          vec augmented = pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data(span::all, span(index));
          bool check = all(to_compare(comparison_inds) == augmented(comparison_inds));
          pvec[static_particle_i].particle_filters[latent_particle_j].consistent(index) = check;
          if(!check) {
            pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data.col(index) =
              initialize_missing_ranks_vec(to_compare, dat.missing_indicator.col(index));
          }
        }

        if(dat.num_new_obs > 0) {
          mat tmp = initialize_missing_ranks(
            dat.new_rankings,
            dat.missing_indicator(
              span::all,
              span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)));

          pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data.resize(dat.n_items, dat.rankings.n_cols);
          pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data(
              span::all,
              span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)
          ) = tmp;
          pvec[static_particle_i].particle_filters[latent_particle_j].log_aug_prob.resize(dat.rankings.n_cols);
        }
      } else if (dat.augpair) {
        for(auto index : dat.updated_match) {
          pvec[static_particle_i].particle_filters[latent_particle_j].consistent(index) = 1;
          vec augmented = pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data(span::all, span(index));
          auto items_above_for_user = dat.items_above[index];
          size_t item1{};
          while(pvec[static_particle_i].particle_filters[latent_particle_j].consistent(index) == 1 && item1 < items_above_for_user.size()) {
            size_t item2{};
            while(pvec[static_particle_i].particle_filters[latent_particle_j].consistent(index) == 1 && item2 < items_above_for_user[item1].size()) {
              if(augmented(items_above_for_user[item1][item2] - 1) > augmented(item1)) {
                pvec[static_particle_i].particle_filters[latent_particle_j].consistent(index) = 0;
              }
              item2++;
            }
            item1++;
          }
        }

        pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data.resize(dat.n_items, dat.rankings.n_cols);

        for(int k{}; k < dat.n_assessors; k++) {
          Rcpp::IntegerVector test = Rcpp::intersect(dat.updated_match, Rcpp::IntegerVector{k});
          if(k >= (dat.n_assessors - dat.num_new_obs) || test.size() > 0) {
            Rcpp::IntegerVector v = Rcpp::sample(sorts[k].n_rows, 1, false, R_NilValue, false);
            ivec ordering = sorts[k].row(v[0]).t();
            uvec rank = sort_index(ordering) + 1;
            pvec[static_particle_i].particle_filters[latent_particle_j].augmented_data.col(k) = conv_to<vec>::from(rank);
          }
        }
        pvec[static_particle_i].particle_filters[latent_particle_j].log_aug_prob.resize(dat.rankings.n_cols);
      }
    }
  }

  return pvec;
}

cube wrapup_rho(const std::vector<std::vector<StaticParticle>>& particle_vectors) {
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

mat wrapup_alpha(const std::vector<std::vector<StaticParticle>>& particle_vectors) {
  mat alpha_samples(particle_vectors[0].size(), particle_vectors.size());
  for(size_t i{}; i < particle_vectors.size(); i++) {
    for(size_t j{}; j < particle_vectors[i].size(); j++) {
      alpha_samples(j, i) = particle_vectors[i][j].alpha;
    }
  }

  return alpha_samples;
}

cube wrapup_augmented_data(const std::vector<StaticParticle>& pvec) {
  cube augmented_data;
  if(!pvec[0].particle_filters[0].augmented_data.is_empty()) {
    augmented_data.set_size(pvec[0].particle_filters[0].augmented_data.n_rows,
                            pvec[0].particle_filters[0].augmented_data.n_cols,
                            pvec.size());
    for(size_t i{}; i < pvec.size(); i++) {
      augmented_data.slice(i) = pvec[i].particle_filters[0].augmented_data;
    }
  }
  return augmented_data;
}

Rcpp::List compute_particle_acceptance(
    const std::vector<std::vector<StaticParticle>>& particle_vectors, int mcmc_steps) {
  vec alpha_acceptance(particle_vectors.size());
  vec rho_acceptance(particle_vectors.size());
  vec aug_acceptance(particle_vectors.size());
  for(size_t i{}; i < particle_vectors.size(); i++) {
    alpha_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const StaticParticle& p) {
          return accumulator + p.alpha_acceptance;
          });
    alpha_acceptance[i] /= particle_vectors[i].size() * mcmc_steps;
    rho_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const StaticParticle& p) {
          return accumulator + p.rho_acceptance;
        });
    rho_acceptance[i] /= particle_vectors[i].size() * mcmc_steps;
    aug_acceptance[i] =
      std::accumulate(
        particle_vectors[i].begin(), particle_vectors[i].end(), 0.0,
        [](double accumulator, const StaticParticle& p) {
          return accumulator + p.particle_filters[0].aug_acceptance / p.particle_filters[0].aug_count;
        });
    aug_acceptance[i] /= particle_vectors[i].size();
  }
  return Rcpp::List::create(
    Rcpp::Named("alpha_acceptance") = alpha_acceptance,
    Rcpp::Named("rho_acceptance") = rho_acceptance,
    Rcpp::Named("aug_acceptance") = aug_acceptance
  );
}
