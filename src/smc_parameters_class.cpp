#include "smc_classes.h"
#include "proposal_functions.h"
using namespace arma;

SMCParameters::SMCParameters(
  const Rcpp::List& model_options,
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
    const Priors& priors) {

  AlphaRatio test = make_new_alpha(
    p.alpha, p.rho,
    alpha_prop_sd, distfun, pfun,
    dat.rankings, dat.observation_frequency,
    dat.n_items, priors
  );
  if(test.accept) p.alpha = test.proposal;
}

void SMCParameters::update_rho(
    Particle& p,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) {
  p.rho = make_new_rho(
    p.rho, dat.rankings, p.alpha, leap_size, distfun,
    dat.observation_frequency);
}

uvec SMCParameters::draw_resampling_index(
  const std::vector<Particle>& pvec
) {
  Rcpp::IntegerVector inds = Rcpp::seq(0, log_inc_wgt.size() - 1);
  vec norm_wgt = exp(log_inc_wgt - max(log_inc_wgt) -
    log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));
  Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(norm_wgt));
  ivec result = Rcpp::sample(inds, log_inc_wgt.size(), true, probs);
  return conv_to<uvec>::from(result);
}
