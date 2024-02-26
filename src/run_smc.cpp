#include <RcppArmadillo.h>

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

  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(pars.metric);
  auto rho_proposal = choose_rank_proposal(
    pars.rho_proposal_option, pars.leap_size);

  dat.update(new_data[0]);
  auto pvec = initialize_particles(initial_values, pars.n_particles, dat);
  aug.reweight(pvec, dat, pfun, distfun);
  pars.resample(pvec);

  par_for_each(
    pvec.begin(), pvec.end(),
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
    Rcpp::Named("rho_samples") = wrapup_rho(pvec),
    Rcpp::Named("alpha_samples") = wrapup_alpha(pvec),
    Rcpp::Named("augmented_rankings") = wrapup_augmented_data(pvec),
    Rcpp::Named("data") = dat.wrapup()
  );

  return particle_history;
}
