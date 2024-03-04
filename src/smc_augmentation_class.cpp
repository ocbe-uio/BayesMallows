#include <limits>
#include "parallel_utils.h"
#include "smc_classes.h"
#include "missing_data.h"
using namespace arma;

unsigned int read_lag(const Rcpp::List& smc_options) {
  Rcpp::IntegerVector tmp = smc_options["latent_sampling_lag"];
  return Rcpp::IntegerVector::is_na(tmp[0]) ?
  std::numeric_limits<unsigned int>::max() :
    static_cast<unsigned int>(tmp[0]);
}

SMCAugmentation::SMCAugmentation(
  const Rcpp::List& compute_options,
  const Rcpp::List& smc_options
  ) :
  partial_aug_prop {
  choose_partial_proposal(compute_options["aug_method"],
                          compute_options["pseudo_aug_metric"]) },
  pairwise_aug_prop {
     choose_pairwise_proposal("none", compute_options["swap_leap"]) },
  latent_sampling_lag { read_lag(smc_options) } {}

void SMCAugmentation::reweight(
    std::vector<Particle>& pvec,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun
) const {
  if(dat.any_missing) {
    par_for_each(
      pvec.begin(), pvec.end(), [distfun = &distfun](Particle& p){
          p.previous_distance = distfun->get()->matdist(p.augmented_data, p.rho);
      });
    augment_partial(pvec, dat);
  } else if(dat.augpair) {
    par_for_each(
      pvec.begin(), pvec.end(), [distfun = &distfun](Particle& p){
        p.previous_distance = distfun->get()->matdist(p.augmented_data, p.rho);
      });
    augment_pairwise(pvec, dat);
  }

  par_for_each(
    pvec.begin(), pvec.end(),
    [n_assessors = dat.n_assessors, num_new_obs = dat.num_new_obs,
     any_missing = dat.any_missing, distfun = &distfun, pfun = &pfun,
     nr = &dat.new_rankings]
    (Particle& p){
      double item_correction_contribution{};
      if(!p.consistent.is_empty()) {
        for(size_t user{}; user < n_assessors - num_new_obs; user++) {
          if(p.consistent(user) == 0) {
            double current_distance = distfun->get()->d(p.augmented_data.col(user), p.rho);

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
          sum(distfun->get()->matdist(new_rankings, p.rho));
      }

      p.log_inc_wgt =
        new_user_contribution + item_correction_contribution -
        num_new_obs * pfun->get()->logz(p.alpha) -
        sum(p.log_aug_prob);
    }
  );
}

void SMCAugmentation::augment_pairwise(
    std::vector<Particle>& pvec, const SMCData& dat
) const {
  par_for_each(
    pvec.begin(), pvec.end(),
    [&dat, pairwise_aug_prop = std::ref(pairwise_aug_prop)]
    (Particle& p){
       for (size_t user{}; user < dat.n_assessors; user++) {
         auto tmp = p.augmented_data.col(user);
         RankProposal pprop = pairwise_aug_prop.get()->propose(
           p.augmented_data.col(user), dat.items_above[user], dat.items_below[user]);
         p.augmented_data.col(user) = pprop.rankings;
         p.log_aug_prob(user) = log(pprop.prob_forward);
       }
     }
  );
}

void SMCAugmentation::augment_partial(
    std::vector<Particle>& pvec, const SMCData& dat
) const {
  par_for_each(
    pvec.begin(), pvec.end(),
    [&dat, partial_aug_prop = std::ref(partial_aug_prop)]
    (Particle& p){
       for (size_t user{}; user < dat.n_assessors; user++) {
        if(user < dat.n_assessors - dat.num_new_obs) {
          if(p.consistent.is_empty()) continue;
          if(p.consistent(user) == 1) continue;
        }

        RankProposal pprop = partial_aug_prop.get()->propose(
          p.augmented_data.col(user), dat.missing_indicator.col(user),
          p.alpha, p.rho);

        p.augmented_data.col(user) = pprop.rankings;
        p.log_aug_prob(user) = log(pprop.prob_forward);
      }
    }
  );
}

void SMCAugmentation::update_missing_ranks(
    Particle& p, const SMCData& dat, const std::unique_ptr<Distance>& distfun) const {
  if(!dat.any_missing && !dat.augpair) return;

  uvec indices_to_loop = find(max(dat.timepoint) - dat.timepoint < latent_sampling_lag);
  for (auto jj : indices_to_loop) {
    std::pair<vec, bool> aug{};
    if(dat.any_missing) {
      aug = make_new_augmentation(
          p.augmented_data.col(jj), dat.missing_indicator.col(jj), p.alpha,
          p.rho, distfun, partial_aug_prop);
    } else if(dat.augpair) {
      aug = make_new_augmentation(
        p.augmented_data.col(jj), p.alpha, p.rho, 0, distfun, pairwise_aug_prop,
        dat.items_above[jj], dat.items_below[jj], "none"
      );
    }
    p.aug_count++;
    if(aug.second) {
      p.augmented_data.col(jj) = aug.first;
      p.aug_acceptance++;
    }

  }
}
