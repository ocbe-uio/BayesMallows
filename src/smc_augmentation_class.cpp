#include "smc_classes.h"
#include "missing_data.h"
using namespace arma;

SMCAugmentation::SMCAugmentation(
  SMCData& dat,
  const Rcpp::List& smc_options,
  const Rcpp::List& initial_values,
  const unsigned int n_particles) :
  aug_method(smc_options["aug_method"]),
  aug_init(initial_values["aug_init"]),
  missing_indicator { set_up_missing(dat) },
  log_aug_prob { arma::zeros(dat.n_assessors, n_particles) }
  {
    if(dat.any_missing){
      augmented_data.set_size(dat.n_items, dat.n_assessors, n_particles);

      for(size_t i{}; i < n_particles; i++) {
        augmented_data.slice(i) =
          initialize_missing_ranks(dat.rankings, missing_indicator);
      }
      if(aug_init.isNotNull()) {
        augmented_data(
          span::all, span(0, dat.rankings.n_cols - dat.num_new_obs - 1), span::all) =
            Rcpp::as<cube>(aug_init);
      }
    }
  }

void SMCAugmentation::reweight(
    SMCParameters& pars,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun
) {
  cube previous_augmented_data = augmented_data;
  augment_partial(pars, dat, distfun);

  for (size_t particle{}; particle < pars.n_particles; ++particle) {
    double item_correction_contribution{};
    if(!dat.consistent.is_empty()) {
      for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
        if(dat.consistent(user, particle) == 0) {
          const arma::vec& pad = previous_augmented_data(span::all, span(user), span(particle));
          const arma::vec& cad = augmented_data(span::all, span(user), span(particle));
          double previous_distance =
            distfun->d(pad, pars.rho_samples.col(particle));
          double current_distance = distfun->d(cad, pars.rho_samples.col(particle));

          item_correction_contribution -=
            pars.alpha_samples(particle) / dat.n_items *
            (current_distance - previous_distance);
        }
      }
    }

    double new_user_contribution{};
    if(dat.num_new_obs > 0) {
      const mat new_rankings = !dat.any_missing ? dat.new_rankings :
      augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1),
        span(particle));
      new_user_contribution = -pars.alpha_samples(particle) / dat.n_items *
        sum(distfun->d(new_rankings, pars.rho_samples.col(particle)));
    }

    pars.log_inc_wgt(particle) =
      new_user_contribution + item_correction_contribution -
      dat.num_new_obs * pfun->logz(pars.alpha_samples(particle)) -
      sum(log_aug_prob.col(particle));
  }
}

void SMCAugmentation::augment_partial(
    const SMCParameters& pars,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun
){
  if(!dat.any_missing) return;
  for (size_t particle{}; particle < pars.n_particles; particle++) {
    for (size_t user{}; user < dat.n_assessors; user++) {
      if(user < dat.n_assessors - dat.num_new_obs) {
        if(dat.consistent.is_empty()) continue;
        if(dat.consistent(user, particle) == 1) continue;
      }

      if (aug_method != "pseudo") {
        augmented_data(span::all, span(user), span(particle)) =
          make_uniform_proposal(
            augmented_data(span::all, span(user), span(particle)),
            missing_indicator.col(user)).rankings;
      } else {
        RankProposal pprop = make_pseudo_proposal(
          augmented_data(span::all, span(user), span(particle)),
          missing_indicator.col(user),
          pars.alpha_samples(particle), pars.rho_samples.col(particle), distfun
        );
        augmented_data(span::all, span(user), span(particle)) = pprop.rankings;
        log_aug_prob(user, particle) = log(pprop.probability);
      }
    }
  }
}

void SMCAugmentation::update_data(
    const unsigned int particle_index, SMCData& dat) {
  if(!dat.any_missing) return;
  dat.rankings = augmented_data.slice(particle_index);
}

void SMCAugmentation::update_missing_ranks(
    const unsigned int particle_index,
    const SMCData& dat,
    const SMCParameters& pars,
    const std::unique_ptr<Distance>& distfun) {
  if(!dat.any_missing) return;

  for (unsigned int jj = dat.n_assessors - dat.num_new_obs;
       jj < dat.n_assessors; ++jj) {
    augmented_data(span::all, span(jj), span(particle_index)) =
      make_new_augmentation(
        augmented_data(span::all, span(jj), span(particle_index)),
        missing_indicator.col(jj),
        pars.alpha_samples(particle_index),
        pars.rho_samples.col(particle_index),
        distfun, log_aug_prob(jj, particle_index), aug_method == "pseudo"
      );
  }
}

void SMCAugmentation::resample(const uvec& index, const SMCData& dat) {
  if(!dat.any_missing) return;
  augmented_data = augmented_data.slices(index);
}
