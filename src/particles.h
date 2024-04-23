#pragma once

#include "smc_classes.h"

std::vector<StaticParticle> initialize_particles(
    const Rcpp::List& initial_values,
    unsigned int n_particles,
    unsigned int n_particle_filters,
    const SMCData& dat
);

std::vector<StaticParticle> prepare_particle_filter(
  const std::vector<StaticParticle>& pvec,
  const SMCData& dat,
  const SMCAugmentation& aug
);

template<typename T>
void resample(std::vector<T>& pvec, const arma::vec& probs,
              const std::unique_ptr<Resampler>& resampler) {
  arma::ivec index = resampler->resample(probs);
  std::vector<T> pvec_old = pvec;
  for(size_t i{}; i < pvec.size(); i++) pvec[i] = pvec_old[index[i]];
  par_for_each(pvec.begin(), pvec.end(), [](T& p) { p.log_inc_wgt = 0;});
}

template<typename T>
arma::vec normalize_probs(std::vector<T>& pvec) {
  arma::vec log_inc_wgt(pvec.size());
  std::transform(pvec.cbegin(), pvec.cend(), log_inc_wgt.begin(),
                 [](const T& p){ return p.log_inc_wgt; });

  return exp(log_inc_wgt - max(log_inc_wgt) -
             log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));
}



arma::mat wrapup_rho(const std::vector<StaticParticle>& particle_vector);
arma::vec wrapup_alpha(const std::vector<StaticParticle>& pvec);
arma::cube wrapup_augmented_data(const std::vector<StaticParticle>& pvec);
