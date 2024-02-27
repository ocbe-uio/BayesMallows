#pragma once

#include "smc_classes.h"

std::vector<Particle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    const SMCData& dat
);

std::vector<Particle> augment_particles(
  const std::vector<Particle>& pvec,
  const SMCData& dat
);

arma::mat wrapup_rho(const std::vector<Particle>& pvec);
arma::vec wrapup_alpha(const std::vector<Particle>& pvec);
arma::cube wrapup_augmented_data(const std::vector<Particle>& pvec);
