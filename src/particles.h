#pragma once

#include "smc_classes.h"

std::vector<Particle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    const SMCData& dat
);

std::vector<Particle> augment_particles(
  const std::vector<Particle>& pvec,
  const SMCData& dat,
  const SMCAugmentation& aug
);

arma::cube wrapup_rho(const std::vector<std::vector<Particle>>& particle_vectors);
arma::mat wrapup_alpha(const std::vector<std::vector<Particle>>& pvec);
arma::cube wrapup_augmented_data(const std::vector<Particle>& pvec);

Rcpp::List compute_particle_acceptance(const std::vector<std::vector<Particle>>& particle_vectors, int mcmc_steps);
