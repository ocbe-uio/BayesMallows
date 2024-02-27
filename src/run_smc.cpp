#include <RcppArmadillo.h>
#include <vector>

#include "parallel_utils.h"
#include "smc_classes.h"
#include "particles.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List  run_smc(
  Rcpp::List data,
  Rcpp::List new_data,
  Rcpp::List model_options,
  Rcpp::List smc_options,
  Rcpp::List compute_options,
  Rcpp::List priors,
  Rcpp::List initial_values,
  Rcpp::Nullable<arma::mat> pfun_values,
  Rcpp::Nullable<arma::mat> pfun_estimate) {

  SMCData dat{data};
  SMCParameters pars{model_options, compute_options, smc_options};
  Priors pris{priors};
  SMCAugmentation aug{compute_options, smc_options};

  std::vector<std::vector<Particle>> particle_vectors(new_data.size() + 1);
  particle_vectors[0] = initialize_particles(initial_values, pars.n_particles, dat);

  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(pars.metric);
  auto rho_proposal = choose_rank_proposal(
    pars.rho_proposal_option, pars.leap_size);



  dat.update(new_data[0]);
  particle_vectors[1] = augment_particles(particle_vectors[0], dat);
  aug.reweight(particle_vectors[1], dat, pfun, distfun);
  pars.resample(particle_vectors[1]);

  par_for_each(
    particle_vectors[1].begin(), particle_vectors[1].end(),
    [&pars, &dat, &pris, &aug, distfun = std::ref(distfun),
     pfun = std::ref(pfun), rho_proposal = std::ref(rho_proposal)]
    (Particle& p) {
       for(size_t i{}; i < pars.mcmc_steps; i++) {
         pars.update_rho(p, dat, distfun, rho_proposal);
         pars.update_alpha(p, dat, pfun, distfun, pris);
         aug.update_missing_ranks(p, dat, distfun);
       }
     }
  );

  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = wrapup_rho(particle_vectors[1]),
    Rcpp::Named("alpha_samples") = wrapup_alpha(particle_vectors[1]),
    Rcpp::Named("augmented_rankings") = wrapup_augmented_data(particle_vectors[1]),
    Rcpp::Named("data") = dat.wrapup()
  );

  return particle_history;
}
