#include <vector>
#include <numeric>
#include "smc_classes.h"
#include "missing_data.h"

using namespace arma;

LatentParticle::LatentParticle(
  const mat& augmented_data, const uvec& particle_consistent,
  const unsigned int n_assessors) :
  augmented_data (augmented_data),
  consistent(particle_consistent),
  log_proposal_prob (zeros(n_assessors)) {}

StaticParticle::StaticParticle(
  double alpha, const vec& rho, const unsigned int n_assessors,
  const std::vector<LatentParticle>& lp) :
  particle_filters (lp),
  alpha (alpha), rho (rho),
  previous_distance(zeros(n_assessors)) {}

void StaticParticle::prepare_particle_filter(const SMCData& dat) {
  alpha_acceptance = 0;
  rho_acceptance = 0;

  for(auto& pf : particle_filters) {
    if(dat.any_missing || dat.augpair) {
      pf.consistent = ones<uvec>(dat.n_assessors - dat.num_new_obs);
    }

    if(dat.any_missing) {
      for(auto index : dat.updated_match) {
        vec to_compare = dat.rankings.col(index);
        uvec comparison_inds = find(to_compare > 0);
        vec augmented = pf.augmented_data(span::all, span(index));
        bool check = all(to_compare(comparison_inds) == augmented(comparison_inds));
        pf.consistent(index) = check;
        if(!check) {
          pf.augmented_data.col(index) =
            initialize_missing_ranks_vec(to_compare, dat.missing_indicator.col(index));
        }
      }

      if(dat.num_new_obs > 0) {
        mat tmp = initialize_missing_ranks(
          dat.new_rankings,
          dat.missing_indicator(
            span::all,
            span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)));

        pf.augmented_data.resize(dat.n_items, dat.rankings.n_cols);
        pf.augmented_data(
            span::all,
            span(dat.rankings.n_cols - dat.num_new_obs, dat.rankings.n_cols - 1)
        ) = tmp;
        pf.log_proposal_prob.resize(dat.rankings.n_cols);
      }
    } else if (dat.augpair) {
      for(auto index : dat.updated_match) {
        pf.consistent(index) = 1;
        vec augmented = pf.augmented_data(span::all, span(index));
        auto items_above_for_user = dat.items_above[index];
        size_t item1{};
        while(pf.consistent(index) == 1 && item1 < items_above_for_user.size()) {
          size_t item2{};
          while(pf.consistent(index) == 1 && item2 < items_above_for_user[item1].size()) {
            if(augmented(items_above_for_user[item1][item2] - 1) > augmented(item1)) {
              pf.consistent(index) = 0;
            }
            item2++;
          }
          item1++;
        }
      }

      pf.augmented_data.resize(dat.n_items, dat.rankings.n_cols);
      pf.log_proposal_prob.resize(dat.rankings.n_cols);
    }
  }
}

std::vector<StaticParticle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    unsigned int n_particle_filters,
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
    latvec.reserve(n_particle_filters);

    for(size_t j{}; j < n_particle_filters; j++) {
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
    }

    pvec.emplace_back(
      StaticParticle(
        alpha_samples(i), rho_samples.col(i), dat.n_assessors,
        latvec));
  }

  return pvec;
}

mat wrapup_rho(const std::vector<StaticParticle>& particle_vector) {
  mat rho_samples(particle_vector[0].rho.size(), particle_vector.size());
  for(size_t j{}; j < particle_vector.size(); j++) {
    rho_samples(span::all, span(j)) = particle_vector[j].rho;
  }
  return rho_samples;
}

vec wrapup_alpha(const std::vector<StaticParticle>& particle_vector) {
  vec alpha_samples(particle_vector.size());
  std::transform(particle_vector.cbegin(), particle_vector.cend(), alpha_samples.begin(),
                [](const StaticParticle& p) { return p.alpha; });

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
