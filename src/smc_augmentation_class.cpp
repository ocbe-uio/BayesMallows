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
    const std::unique_ptr<Distance>& distfun
) {
  if(dat.any_missing) {
    std::for_each(
      pvec.begin(), pvec.end(),
      [distfun = &distfun](Particle& p){
          p.previous_distance = distfun->get()->d(p.augmented_data, p.rho);
      }
    );
  }

  if(dat.any_missing) augment_partial(pvec, dat);

  for (size_t particle{}; particle < pvec.size(); ++particle) {
    double item_correction_contribution{};
    if(!pvec[particle].consistent.is_empty()) {
      for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
        if(pvec[particle].consistent(user) == 0) {
          const arma::vec& cad = pvec[particle].augmented_data.col(user);
          double current_distance = distfun->d(cad, pvec[particle].rho);

          item_correction_contribution -= pvec[particle].alpha / dat.n_items *
            (current_distance - pvec[particle].previous_distance(user));
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
    const SMCData& dat
){
  std::for_each(
    pvec.begin(), pvec.end(),
    [n_assessors = dat.n_assessors, num_new_obs = dat.num_new_obs,
     aug_method = aug_method, pseudo_aug_metric = pseudo_aug_metric,
     missing_indicator = missing_indicator]
    (Particle& p){
       auto pseudo_aug_distance = aug_method == "uniform" ? nullptr :
       choose_distance_function(pseudo_aug_metric);

       for (size_t user{}; user < n_assessors; user++) {
        if(user < n_assessors - num_new_obs) {
          if(p.consistent.is_empty()) continue;
          if(p.consistent(user) == 1) continue;
        }
        if (pseudo_aug_distance == nullptr) {
          p.augmented_data.col(user) =
            make_uniform_proposal(
              p.augmented_data.col(user),
              missing_indicator.col(user)).rankings;
        } else {
          RankProposal pprop = make_pseudo_proposal(
            p.augmented_data.col(user),
            missing_indicator.col(user),
            p.alpha, p.rho,
            pseudo_aug_distance
          );
          p.augmented_data.col(user) = pprop.rankings;
          p.log_aug_prob(user) = log(pprop.probability);
        }
      }
    }
  );
}

void SMCAugmentation::update_missing_ranks(
    Particle& p,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) {
  if(!dat.any_missing) return;

  auto pseudo_aug_distance = aug_method == "uniform" ? nullptr : choose_distance_function(pseudo_aug_metric);

  for (unsigned int jj{}; jj < dat.n_assessors; ++jj) {
    p.augmented_data.col(jj) =
      make_new_augmentation(
        p.augmented_data.col(jj), missing_indicator.col(jj), p.alpha,
        p.rho, distfun, pseudo_aug_distance, p.log_aug_prob(jj));
  }
}
