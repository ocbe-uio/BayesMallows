#include "classes.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "distances.h"
#include "partitionfuns.h"
using namespace arma;

Augmentation::Augmentation(
  Data& dat,
  const Rcpp::List& compute_options
) :
  augpair { dat.constraints.length() > 0 },
  any_missing { !is_finite(dat.rankings) },
  save_aug { compute_options["save_aug"] },
  aug_thinning { compute_options["aug_thinning"] },
  swap_leap { compute_options["swap_leap"] } {

    if(any_missing){
      set_up_missing(dat.rankings, missing_indicator);
      initialize_missing_ranks(dat.rankings, missing_indicator);
    }
    if(save_aug){
      unsigned int nmc = Rcpp::as<unsigned int>(compute_options["nmc"]);
      augmented_data.set_size(dat.n_items, dat.n_assessors,
                              std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
      augmented_data.slice(0) = dat.rankings;
    }}

SMCAugmentation::SMCAugmentation(
  SMCData& dat,
  const Rcpp::List& smc_options,
  const Rcpp::List& initial_values,
  const unsigned int n_particles) :
  aug_method(smc_options["aug_method"]),
  log_aug_prob { arma::zeros(n_particles) },
  any_missing { !is_finite(dat.rankings) },
  aug_init(initial_values["aug_init"])
  {
    if(any_missing){
      set_up_missing(dat.rankings, missing_indicator);
      augmented_data.set_size(dat.n_items, dat.n_assessors, n_particles);

      for(size_t i{}; i < n_particles; i++) {
        augmented_data.slice(i) = dat.rankings;
        initialize_missing_ranks(augmented_data.slice(i), missing_indicator);
      }
      if(aug_init.isNotNull()) {
        augmented_data(
          span::all,
          span(0, dat.rankings.n_cols - dat.new_rankings.n_cols - 1),
          span::all) = Rcpp::as<cube>(aug_init);
      }
    }
  }

void Augmentation::augment_pairwise(
    const unsigned int t,
    Data& dat,
    const Parameters& pars,
    const Clustering& clus
){
  if(!augpair) return;
  for(size_t i = 0; i < dat.n_assessors; ++i) {
    vec proposal;
    int g_diff = 0;
    if(pars.error_model == "none"){
      proposal = propose_pairwise_augmentation(dat.rankings.col(i), Rcpp::as<Rcpp::List>(dat.constraints[i]));
    } else if(pars.error_model == "bernoulli"){
      proposal = propose_swap(dat.rankings.col(i), Rcpp::as<Rcpp::List>(dat.constraints[i]), g_diff, swap_leap);
    } else {
      Rcpp::stop("error_model must be 'none' or 'bernoulli'");
    }

    double u = std::log(randu<double>());

    // Find which cluster the assessor belongs to
    int cluster = clus.current_cluster_assignment(i);

    double ratio = -pars.alpha_old(cluster) / dat.n_items *
      (get_rank_distance(proposal, pars.rho_old.col(cluster), pars.metric) -
      get_rank_distance(dat.rankings.col(i), pars.rho_old.col(cluster), pars.metric));

    if(pars.error_model != "none") {
      ratio += g_diff * std::log(pars.theta(t) / (1 - pars.theta(t)));
    }

    if(ratio > u){
      dat.rankings.col(i) = proposal;
    }
  }
}

void Augmentation::update_missing_ranks(
    Data& dat,
    const Clustering& clus,
    const Parameters& pars) {
  if(!any_missing) return;

  for(size_t i = 0; i < dat.n_assessors; ++i){
    int cluster = clus.current_cluster_assignment(i);

    dat.rankings.col(i) = make_new_augmentation(
      dat.rankings.col(i), missing_indicator.col(i),
      pars.alpha_old(cluster), pars.rho_old.col(cluster),
      pars.metric
    );
  }
}

void SMCAugmentation::reweight(
    SMCParameters& pars,
    const SMCData& dat,
    const Rcpp::List& logz_list
) {
  augment_partial(pars, dat);

  for (size_t particle{}; particle < pars.n_particles; ++particle) {
    const double log_z_alpha = get_partition_function(
      dat.n_items, pars.alpha_samples(particle), logz_list, pars.metric
    );
    double item_correction_contribution{};
    if(!dat.consistent.is_empty()) {
      for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
        if(dat.consistent(user, particle) == 0) {
          vec previous_augmented_ranking = augmented_data(span::all, span(user), span(particle - 1));
          vec current_augmented_ranking = augmented_data(span::all, span(user), span(particle));

          double previous_distance = get_rank_distance(
            previous_augmented_ranking, pars.rho_samples.col(particle),
            pars.metric);
          double current_distance = get_rank_distance(
            current_augmented_ranking, pars.rho_samples.col(particle),
            pars.metric);

          item_correction_contribution -= pars.alpha_samples(particle) / dat.n_items *
            (current_distance - previous_distance);
        }
      }
    }

    double new_user_contribution{};
    if(dat.num_new_obs > 0) {
      const mat new_rankings = !any_missing ? dat.new_rankings :
      augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1),
        span(particle));
      new_user_contribution = -pars.alpha_samples(particle) / dat.n_items *
        rank_dist_sum(new_rankings, pars.rho_samples.col(particle), pars.metric,
                      dat.observation_frequency(span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1)));
    }

    pars.log_inc_wgt(particle) =
      new_user_contribution + item_correction_contribution -
      dat.num_new_obs * log_z_alpha - log_aug_prob(particle);
  }
}

void SMCAugmentation::augment_partial(
    const SMCParameters& pars,
    const SMCData& dat
){
  if(!any_missing) return;
  for (size_t particle{}; particle < pars.n_particles; ++particle) {

    for (size_t user{}; user < dat.n_assessors; ++user) {
      if(user < dat.n_assessors - dat.num_new_obs) {
        if(dat.consistent.is_empty()) continue;
        if(dat.consistent(user, particle) == 1) continue;
      }

      uvec unranked_items = shuffle(find(missing_indicator.col(user) == 1));

      if (aug_method != "pseudo") {
        augmented_data(span::all, span(user), span(particle)) =
          propose_augmentation(augmented_data(span::all, span(user), span(particle)),
                               missing_indicator.col(user));

      } else {
        PseudoProposal pprop = make_pseudo_proposal(
          unranked_items, augmented_data(span::all, span(user), span(particle)),
          pars.alpha_samples(particle), pars.rho_samples.col(particle), pars.metric
        );
        augmented_data(span::all, span(user), span(particle)) = pprop.rankings;
        log_aug_prob(particle) += log(pprop.probability);
      }
    }
  }
}

void SMCAugmentation::update_data(
    const unsigned int particle_index, SMCData& dat) {
  if(!any_missing) return;
  dat.rankings = augmented_data.slice(particle_index);
}
