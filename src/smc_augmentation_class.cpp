#include <limits>
#include <numeric>
#include "parallel_utils.h"
#include "particles.h"
#include "smc_classes.h"
#include "missing_data.h"
#include "all_topological_sorts.h"
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
  max_topological_sorts{ smc_options["max_topological_sorts"] },
  partial_aug_prop {
  choose_partial_proposal(compute_options["aug_method"],
                          compute_options["pseudo_aug_metric"]) },
  pairwise_aug_prop {
     choose_pairwise_proposal("none", compute_options["swap_leap"]) },
  latent_sampling_lag { read_lag(smc_options) } {}

void SMCAugmentation::run_particle_filter(
    std::vector<StaticParticle>& pvec,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<Resampler>& resampler,
    size_t time
) const {

  for(auto& p : pvec){
    if(time > 0) {
      vec resampling_probs = normalize_probs(p.particle_filters);
      resample(p.particle_filters, resampling_probs, resampler);
    }
    for(auto& pf : p.particle_filters) {
      if(dat.augpair) {
        for(int k{}; k < dat.n_assessors; k++) {
          Rcpp::IntegerVector test = Rcpp::intersect(dat.updated_match, Rcpp::IntegerVector{k});
          if(k >= (dat.n_assessors - dat.num_new_obs) || test.size() > 0) {
            uvec indices = find(dat.preferences.col(0) == dat.user_ids[k]);
            imat prefs = dat.preferences.rows(indices);
            imat sorts = all_topological_sorts(prefs.cols(1, 2), dat.n_items, 1e5);

            Rcpp::IntegerVector v = Rcpp::sample(sorts.n_rows, 1, false, R_NilValue, false);
            ivec ordering = sorts.row(v[0]).t();
            uvec rank = sort_index(ordering) + 1;
            pf.augmented_data.col(k) = conv_to<vec>::from(rank);
          }
        }
      }

      if(dat.any_missing || dat.augpair) {
        p.previous_distance = distfun->matdist(pf.augmented_data, p.rho);
        p = augment_partial(p, dat);
      }

      double item_correction_contribution{};
      if(!pf.consistent.is_empty()) {
        for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
          if(pf.consistent(user) == 0) {
            double current_distance = distfun->d(pf.augmented_data.col(user), p.rho);

            item_correction_contribution -= p.alpha / p.rho.size() *
              (current_distance - p.previous_distance(user));
          }
        }
      }

      double new_user_contribution{};
      if(dat.num_new_obs > 0) {
        mat new_rankings;
        if(dat.any_missing || dat.augpair) {
          new_rankings = pf.augmented_data(
            span::all,
            span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1));
        } else {
          new_rankings = dat.new_rankings;
        }

        new_user_contribution = -p.alpha / p.rho.size() *
          sum(distfun->matdist(new_rankings, p.rho));
      }

      pf.log_inc_wgt =
        new_user_contribution + item_correction_contribution -
        dat.num_new_obs * pfun->logz(p.alpha) -
        sum(pf.log_proposal_prob);
    }

    vec log_weights(p.particle_filters.size());
    std::transform(
      p.particle_filters.cbegin(), p.particle_filters.cend(), log_weights.begin(),
      [](const LatentParticle& lp){ return lp.log_inc_wgt; });

    double log_average_weights = max(log_weights) + log(sum(exp(log_weights - max(log_weights)))) -
      log(p.particle_filters.size());

    // log_inc_wgt is reset under resampling but marginal_log_likelihood is not.
    p.marginal_log_likelihood += log_average_weights;
    p.log_inc_wgt += log_average_weights;
  }
}

StaticParticle SMCAugmentation::augment_partial(
    const StaticParticle& p, const SMCData& dat
) const {
  StaticParticle ret{p};

  for (size_t user{}; user < dat.n_assessors; user++) {
  if(user < dat.n_assessors - dat.num_new_obs) {
    if(ret.particle_filters[0].consistent.is_empty()) continue;
    if(ret.particle_filters[0].consistent(user) == 1) continue;
  }

  RankProposal pprop;
  if(dat.any_missing) {
    pprop = partial_aug_prop.get()->propose(
      ret.particle_filters[0].augmented_data.col(user), dat.missing_indicator.col(user),
      ret.alpha, p.rho);
  } else if(dat.augpair) {
    pprop = pairwise_aug_prop.get()->propose(
      ret.particle_filters[0].augmented_data.col(user), dat.items_above[user], dat.items_below[user]);
  }

  ret.particle_filters[0].augmented_data.col(user) = pprop.rankings;
  ret.particle_filters[0].log_proposal_prob(user) = log(pprop.prob_forward);
  }
  return ret;
}

void SMCAugmentation::update_missing_ranks(
    StaticParticle& p, const SMCData& dat, const std::unique_ptr<Distance>& distfun) const {
  if(!dat.any_missing && !dat.augpair) return;

  uvec indices_to_loop = find(max(dat.timepoint) - dat.timepoint < latent_sampling_lag);
  for (auto jj : indices_to_loop) {
    std::pair<vec, bool> aug{};
    if(dat.any_missing) {
      aug = make_new_augmentation(
          p.particle_filters[0].augmented_data.col(jj), dat.missing_indicator.col(jj), p.alpha,
          p.rho, distfun, partial_aug_prop);
    } else if(dat.augpair) {
      aug = make_new_augmentation(
        p.particle_filters[0].augmented_data.col(jj), p.alpha, p.rho, 0, distfun, pairwise_aug_prop,
        dat.items_above[jj], dat.items_below[jj], "none"
      );
    }
    p.particle_filters[0].aug_count++;
    if(aug.second) {
      p.particle_filters[0].augmented_data.col(jj) = aug.first;
      p.particle_filters[0].aug_acceptance++;
    }

  }
}
