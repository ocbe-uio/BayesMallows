#pragma once

#include "smc_classes.h"

std::vector<StaticParticle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    const SMCData& dat
);

std::vector<StaticParticle> augment_particles(
  const std::vector<StaticParticle>& pvec,
  const SMCData& dat,
  const SMCAugmentation& aug
);

arma::cube wrapup_rho(const std::vector<std::vector<StaticParticle>>& particle_vectors);
arma::mat wrapup_alpha(const std::vector<std::vector<StaticParticle>>& pvec);
arma::cube wrapup_augmented_data(const std::vector<StaticParticle>& pvec);

Rcpp::List compute_particle_acceptance(const std::vector<std::vector<StaticParticle>>& particle_vectors, int mcmc_steps);
