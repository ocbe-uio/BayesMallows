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

  SMCData dat{data, new_data};
  SMCParameters pars{compute_options, smc_options};
  Priors pris{priors};
  SMCAugmentation aug{compute_options, smc_options};
  std::string metric = model_options["metric"];
  auto pfun = choose_partition_function(
    dat.n_items, metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(metric);
  auto pvec = initialize_particles(data, initial_values, smc_options, aug, dat);
  auto rho_proposal = choose_rank_proposal(
    pars.rho_proposal_option, pars.leap_size);
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
    Rcpp::Named("augmented_rankings") = wrapup_augmented_data(pvec)
  );

  return particle_history;
}
