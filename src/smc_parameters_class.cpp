#include <utility>
#include "smc_classes.h"
#include "proposal_functions.h"
#include "parallel_utils.h"
using namespace arma;

SMCParameters::SMCParameters(
  const Rcpp::List& model_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& smc_options
) :
  mcmc_steps { smc_options["mcmc_steps"] },
  metric ( model_options["metric"] ),
  n_particles { smc_options["n_particles"] },
  n_particle_filters { smc_options["n_particle_filters"] },
  resampling_threshold { smc_options["resampling_threshold"]},
  rho_proposal_function {
    choose_rho_proposal(compute_options["rho_proposal"],
                        compute_options["leap_size"])},
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  resampling_method(smc_options["resampler"]) {}

void SMCParameters::update_alpha(
    StaticParticle& p,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const Priors& priors) const {

  AlphaRatio test = make_new_alpha(
    p.alpha, p.rho, alpha_prop_sd, distfun, pfun,
    !p.particle_filters[0].augmented_data.empty() ? p.particle_filters[0].augmented_data : dat.rankings,
    dat.observation_frequency, dat.n_items, priors
  );
  if(test.accept) {
    p.alpha = test.proposal;
    p.alpha_acceptance++;
  }
}

void SMCParameters::update_rho(
    StaticParticle& p,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) const {
  auto proposal = make_new_rho(
    p.rho, !p.particle_filters[0].augmented_data.empty() ? p.particle_filters[0].augmented_data : dat.rankings,
    p.alpha, distfun, rho_proposal_function, dat.observation_frequency);
  if(proposal.second) {
    p.rho = proposal.first;
    p.rho_acceptance++;
  }
}

