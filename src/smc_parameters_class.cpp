#include "smc_classes.h"
#include "proposal_functions.h"
using namespace arma;

SMCParameters::SMCParameters(
  const Rcpp::List& model_options,
  const Rcpp::List& smc_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values
) :
  n_particles { smc_options["n_particles"] },
  mcmc_steps { smc_options["mcmc_steps"] },
  alpha_samples(initial_values["alpha_init"]) ,
  rho_samples(initial_values["rho_init"]),
  log_inc_wgt { zeros(n_particles) },
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  leap_size { compute_options["leap_size"] } {}

void SMCParameters::update_rho(
    const unsigned int particle_index,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) {
  rho_samples.col(particle_index) = make_new_rho(
    rho_samples.col(particle_index), dat.rankings,
    alpha_samples(particle_index), leap_size, distfun,
    dat.observation_frequency);
}

void SMCParameters::update_alpha(
    const unsigned int particle_index,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const Priors& priors) {

  AlphaRatio test = make_new_alpha(
    alpha_samples(particle_index),
    rho_samples.col(particle_index),
    alpha_prop_sd, distfun, pfun,
    dat.rankings, dat.observation_frequency,
    dat.n_items, priors
  );
  if(test.accept){
    alpha_samples(particle_index) = test.proposal;
  }
}

uvec SMCParameters::draw_resampling_index() {
  Rcpp::IntegerVector inds = Rcpp::seq(0, log_inc_wgt.size() - 1);
  vec norm_wgt = exp(log_inc_wgt - max(log_inc_wgt) -
    log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));
  Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(norm_wgt));
  ivec result = Rcpp::sample(inds, log_inc_wgt.size(), true, probs);
  return conv_to<uvec>::from(result);
}
