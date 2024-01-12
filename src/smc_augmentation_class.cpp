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
      pvec.begin(), pvec.end(), [distfun = &distfun](Particle& p){
          p.previous_distance = distfun->get()->d(p.augmented_data, p.rho);
      });
    augment_partial(pvec, dat);
  }

  std::for_each(
    pvec.begin(), pvec.end(),
    [n_assessors = dat.n_assessors, num_new_obs = dat.num_new_obs,
     any_missing = dat.any_missing, distfun = &distfun, pfun = &pfun,
     nr = &dat.new_rankings]
    (Particle& p){
      double item_correction_contribution{};
      if(!p.consistent.is_empty()) {
        for(size_t user{}; user < n_assessors - num_new_obs; user++) {
          if(p.consistent(user) == 0) {
            const arma::vec& cad = p.augmented_data.col(user);
            double current_distance = distfun->get()->d(cad, p.rho);

            item_correction_contribution -= p.alpha / p.rho.size() *
              (current_distance - p.previous_distance(user));
          }
        }
      }

      double new_user_contribution{};
      if(num_new_obs > 0) {
        mat new_rankings;
        if(any_missing) {
          new_rankings = p.augmented_data(
            span::all,
            span(n_assessors - num_new_obs, n_assessors - 1));
        } else {
          new_rankings = *nr;
        }

        new_user_contribution = -p.alpha / p.rho.size() *
          sum(distfun->get()->d(new_rankings, p.rho));
      }

      p.log_inc_wgt =
        new_user_contribution + item_correction_contribution -
        num_new_obs * pfun->get()->logz(p.alpha) -
        sum(p.log_aug_prob);
    }
  );

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
