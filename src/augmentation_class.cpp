#include "classes.h"
#include "missing_data.h"
#include "distances.h"
#include "rank_proposal.h"
using namespace arma;

Augmentation::Augmentation(
  Data& dat,
  const Rcpp::List& compute_options,
  const Rcpp::List& model_options
) :
  save_aug { compute_options["save_aug"] },
  aug_thinning { compute_options["aug_thinning"] },
  error_model(model_options["error_model"]),
  partial_aug_prop {
    choose_partial_proposal(compute_options["aug_method"],
                            compute_options["pseudo_aug_metric"]) },
  pairwise_aug_prop {
    choose_pairwise_proposal(error_model, compute_options["swap_leap"])
  }
   {
    if(dat.any_missing){
      dat.rankings = initialize_missing_ranks(dat.rankings, dat.missing_indicator);
    }
    if(save_aug){
      unsigned int nmc{ compute_options["nmc"] };
      augmented_data.set_size(
        dat.n_items, dat.n_assessors,
        std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
      augmented_data.slice(0) = dat.rankings;
    }}

void Augmentation::update_missing_ranks(
    Data& dat,
    const Clustering& clus,
    const Parameters& pars,
    const std::unique_ptr<Distance>& distfun) {
  if(!dat.any_missing && !dat.augpair) return;
  for(size_t i = 0; i < dat.n_assessors; ++i){
    int cluster = clus.current_cluster_assignment(i);
    std::pair<vec, bool> aug{};
    if(dat.any_missing) {
      aug = make_new_augmentation(
        dat.rankings.col(i), dat.missing_indicator.col(i),
        pars.alpha_old(cluster), pars.rho_old.col(cluster),
        distfun, partial_aug_prop);
    } else if(dat.augpair) {
      aug = make_new_augmentation(
        dat.rankings.col(i), pars.alpha_old(cluster), pars.rho_old.col(cluster),
        pars.theta_current, distfun, pairwise_aug_prop,
        dat.items_above[i], dat.items_below[i], error_model);
    }
    if(pars.t > pars.burnin) aug_count++;
    if(aug.second) {
      dat.rankings.col(i) = aug.first;
      if(pars.t > pars.burnin) aug_acceptance++;
    }
  }
}

void Augmentation::save_augmented_data(const Data& dat, const Parameters& pars) {
  if(save_aug & (pars.t % aug_thinning == 0)){
    ++aug_index;
    augmented_data.slice(aug_index) = dat.rankings;
  }
}
