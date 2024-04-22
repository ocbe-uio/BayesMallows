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
  rho_proposal_function {
    choose_rho_proposal(compute_options["rho_proposal"],
                        compute_options["leap_size"])},
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  resampler { choose_resampler(std::string(smc_options["resampler"])) },
  resampling_threshold { smc_options["resampling_threshold"]} {}

void SMCParameters::update_alpha(
    StaticParticle& p,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const Priors& priors) const {

  AlphaRatio test = make_new_alpha(
    p.alpha, p.rho, alpha_prop_sd, distfun, pfun,
    !p.augmented_data.empty() ? p.augmented_data : dat.rankings,
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
    p.rho, !p.augmented_data.empty() ? p.augmented_data : dat.rankings,
    p.alpha, distfun, rho_proposal_function, dat.observation_frequency);
  if(proposal.second) {
    p.rho = proposal.first;
    p.rho_acceptance++;
  }
}

void SMCParameters::resample(std::vector<StaticParticle>& pvec) {

  arma::vec log_inc_wgt(pvec.size());
  std::transform(pvec.cbegin(), pvec.cend(), log_inc_wgt.begin(),
                 [](const StaticParticle& p){ return p.log_inc_wgt; });

  arma::vec probs = exp(log_inc_wgt - max(log_inc_wgt) -
    log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));

  double ess = 1 / accu(pow(probs, 2));

  if(ess < resampling_threshold) {
    rejuvenate = true;
    par_for_each(pvec.begin(), pvec.end(), [](StaticParticle& p) { p.log_inc_wgt = 0;});
    arma::ivec index = resampler->resample(probs);
    std::vector<StaticParticle> pvec_old = pvec;
    for(size_t i{}; i < pvec.size(); i++) pvec[i] = pvec_old[index[i]];
  } else {
    rejuvenate = false;
  }
}
