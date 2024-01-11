#include "smc_classes.h"

using namespace arma;

Particle::Particle(double alpha, const vec& rho) :
  alpha (alpha), rho (rho) {}

std::vector<Particle> initialize_particles(
    const Rcpp::List& initial_values,
    const Rcpp::List& smc_options
) {
  const unsigned int n_particles { smc_options["n_particles"] };
  vec alpha_samples(initial_values["alpha_init"]);
  mat rho_samples(initial_values["rho_init"]);

  std::vector<Particle> pvec;
  pvec.reserve(n_particles);

  for(size_t i{}; i < n_particles; i++) {
    pvec.emplace_back(
      Particle(alpha_samples(i), rho_samples.col(i))
    );
  }

  return pvec;
}
