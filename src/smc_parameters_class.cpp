#include "smc_classes.h"
#include "proposal_functions.h"
using namespace arma;

SMCParameters::SMCParameters(
  const Rcpp::List& smc_options,
  const Rcpp::List& compute_options
) :
  mcmc_steps { smc_options["mcmc_steps"] },
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  leap_size { compute_options["leap_size"] } {}

void SMCParameters::update_alpha(
    Particle& p,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const Priors& priors) const {

  AlphaRatio test = make_new_alpha(
    p.alpha, p.rho, alpha_prop_sd, distfun, pfun,
    dat.any_missing ? p.augmented_data : dat.rankings,
    dat.observation_frequency, dat.n_items, priors
  );
  if(test.accept) p.alpha = test.proposal;
}

void SMCParameters::update_rho(
    Particle& p,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) const {
  p.rho = make_new_rho(
    p.rho, dat.any_missing ? p.augmented_data : dat.rankings,
    p.alpha, leap_size, distfun,
    dat.observation_frequency);
}

Rcpp::IntegerVector SMCParameters::draw_resampling_index(
  const std::vector<Particle>& pvec
) const {
  Rcpp::NumericVector log_inc_wgt(pvec.size());
  std::transform(pvec.cbegin(), pvec.cend(), log_inc_wgt.begin(),
                 [](const Particle& p){ return p.log_inc_wgt; });

  Rcpp::NumericVector probs = exp(log_inc_wgt - max(log_inc_wgt) -
    log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));

  return Rcpp::sample(log_inc_wgt.size(), log_inc_wgt.size(), true, probs, false);
}

void SMCParameters::resample(std::vector<Particle>& pvec) const {
  Rcpp::IntegerVector index = draw_resampling_index(pvec);
  std::vector<Particle> pvec_old = pvec;
  for(size_t i{}; i < pvec.size(); i++) pvec[i] = pvec_old[index[i]];
}
