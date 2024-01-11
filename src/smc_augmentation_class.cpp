#include "smc_classes.h"
#include "missing_data.h"
using namespace arma;

SMCAugmentation::SMCAugmentation(
  SMCData& dat,
  const Rcpp::List& compute_options) :
  missing_indicator { set_up_missing(dat) },
  aug_method(compute_options["aug_method"]),
  pseudo_aug_metric(compute_options["pseudo_aug_metric"]) {}

void SMCAugmentation::reweight(
    std::vector<Particle>& pvec,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<Distance>& pseudo_aug_distance
) {
  cube previous_augmented_data;
  if(dat.any_missing) {
    previous_augmented_data.set_size(dat.n_items, dat.n_assessors, pvec.size());
    for(size_t i{}; i < pvec.size(); i++) {
      previous_augmented_data.slice(i) = pvec[i].augmented_data;
    }
  }

  if(dat.any_missing) augment_partial(pvec, dat, pseudo_aug_distance);

  for (size_t particle{}; particle < pvec.size(); ++particle) {
    double item_correction_contribution{};
    if(!pvec[particle].consistent.is_empty()) {
      for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
        if(pvec[particle].consistent(user) == 0) {
          const arma::vec& pad = previous_augmented_data(span::all, span(user), span(particle));
          const arma::vec& cad = pvec[particle].augmented_data.col(user);
          double previous_distance =
            distfun->d(pad, pvec[particle].rho);
          double current_distance = distfun->d(cad, pvec[particle].rho);

          item_correction_contribution -= pvec[particle].alpha / dat.n_items *
            (current_distance - previous_distance);
        }
      }
    }

    double new_user_contribution{};
    if(dat.num_new_obs > 0) {
      const mat new_rankings = !dat.any_missing ? dat.new_rankings :
      pvec[particle].augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1));
      new_user_contribution = -pvec[particle].alpha / dat.n_items *
        sum(distfun->d(new_rankings, pvec[particle].rho));
    }

    pvec[particle].log_inc_wgt =
      new_user_contribution + item_correction_contribution -
      dat.num_new_obs * pfun->logz(pvec[particle].alpha) -
      sum(pvec[particle].log_aug_prob);
  }
}

void SMCAugmentation::augment_partial(
    std::vector<Particle>& pvec,
    const SMCData& dat,
    const std::unique_ptr<Distance>& pseudo_aug_distance
){
  // std::for_each(
  //   pvec.begin(), pvec.end(),
  //   [n_assessors = dat.n_assessors, num_new_obs = dat.num_new_obs]
  //   (Particle& p){
  //     for (size_t user{}; user < n_assessors; user++) {
  //       if(user < n_assessors - num_new_obs) {
  //         if(p.consistent.is_empty()) continue;
  //         if(p.consistent(user) == 1) continue;
  //       }
  //       if (pseudo_aug_distance == nullptr) {
  //         p.augmented_data.col(user) =
  //           make_uniform_proposal(
  //             p.augmented_data.col(user),
  //             missing_indicator.col(user)).rankings;
  //       } else {
  //         RankProposal pprop = make_pseudo_proposal(
  //           p.augmented_data.col(user),
  //           missing_indicator.col(user),
  //           p.alpha, p.rho,
  //           pseudo_aug_distance
  //         );
  //         p.augmented_data.col(user) = pprop.rankings;
  //         p.log_aug_prob(user) = log(pprop.probability);
  //       }
  //     }
  //   }
  // )
  for (size_t particle{}; particle < pvec.size(); particle++) {
    for (size_t user{}; user < dat.n_assessors; user++) {
      if(user < dat.n_assessors - dat.num_new_obs) {
        if(pvec[particle].consistent.is_empty()) continue;
        if(pvec[particle].consistent(user) == 1) continue;
      }
      if (pseudo_aug_distance == nullptr) {
        pvec[particle].augmented_data.col(user) =
          make_uniform_proposal(
            pvec[particle].augmented_data.col(user),
            missing_indicator.col(user)).rankings;
      } else {
        RankProposal pprop = make_pseudo_proposal(
          pvec[particle].augmented_data.col(user),
          missing_indicator.col(user),
          pvec[particle].alpha, pvec[particle].rho,
          pseudo_aug_distance
        );
        pvec[particle].augmented_data.col(user) = pprop.rankings;
        pvec[particle].log_aug_prob(user) = log(pprop.probability);
      }
    }
  }
}

void SMCAugmentation::update_missing_ranks(
    Particle& p,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<Distance>& pseudo_aug_distance) {
  if(!dat.any_missing) return;

  for (unsigned int jj{}; jj < dat.n_assessors; ++jj) {
    p.augmented_data.col(jj) =
      make_new_augmentation(
        p.augmented_data.col(jj), missing_indicator.col(jj), p.alpha,
        p.rho, distfun, pseudo_aug_distance, p.log_aug_prob(jj));
  }
}
