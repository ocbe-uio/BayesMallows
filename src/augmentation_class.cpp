#include "classes.h"
#include "missing_data.h"
#include "distances.h"
#include "rank_proposal.h"
using namespace arma;

Augmentation::Augmentation(
  Data& dat,
  const Rcpp::List& compute_options
) :
  augpair { dat.items_above.size() > 0 },
  save_aug { compute_options["save_aug"] },
  aug_thinning { compute_options["aug_thinning"] },
  swap_leap { compute_options["swap_leap"] },
  missing_indicator { set_up_missing(dat) },
  aug_method ( compute_options["aug_method"] ),
  pseudo_aug_metric ( compute_options["pseudo_aug_metric"] ),
  pseudo_aug_distance {
    aug_method == "uniform" ? nullptr : choose_distance_function(pseudo_aug_metric)
  },
  log_aug_prob { zeros(dat.n_assessors) } {
    if(dat.any_missing){
      dat.rankings = initialize_missing_ranks(dat.rankings, missing_indicator);
    }
    if(save_aug){
      unsigned int nmc{ compute_options["nmc"] };
      augmented_data.set_size(
        dat.n_items, dat.n_assessors,
        std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
      augmented_data.slice(0) = dat.rankings;
    }}

void Augmentation::augment_pairwise(
    const unsigned int t,
    Data& dat,
    const Parameters& pars,
    const Clustering& clus,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<ProposalDistribution>& prop
){
  if(!augpair) return;
  for(size_t i = 0; i < dat.n_assessors; ++i) {

    RankProposal rp = prop->propose(
        dat.rankings.col(i), dat.items_above[i], dat.items_below[i]);

    double u = std::log(R::runif(0, 1));
    int cluster = clus.current_cluster_assignment(i);

    double newdist = distfun->d(rp.rankings, pars.rho_old.col(cluster), rp.mutated_items);
    double olddist = distfun->d(dat.rankings.col(i), pars.rho_old.col(cluster), rp.mutated_items);
    double ratio = -pars.alpha_old(cluster) / dat.n_items * (newdist - olddist);

    if(pars.error_model != "none") {
      ratio += rp.g_diff * std::log(pars.theta(t) / (1 - pars.theta(t)));
    }

    if(ratio > u) dat.rankings.col(i) = rp.rankings;
  }
}

void Augmentation::update_missing_ranks(
    Data& dat,
    const Clustering& clus,
    const Parameters& pars,
    const std::unique_ptr<Distance>& distfun) {
  if(!dat.any_missing) return;

  for(size_t i = 0; i < dat.n_assessors; ++i){
    int cluster = clus.current_cluster_assignment(i);
    dat.rankings.col(i) = make_new_augmentation(
      dat.rankings.col(i), missing_indicator.col(i),
      pars.alpha_old(cluster), pars.rho_old.col(cluster),
      distfun, pseudo_aug_distance,
      log_aug_prob(i)
    );
  }
}
