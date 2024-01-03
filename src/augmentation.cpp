#include "classes.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "distances.h"
using namespace arma;

Augmentation::Augmentation(
  Data& dat,
  const Rcpp::List& compute_options
) :
  augpair { dat.items_above.size() > 0 },
  save_aug { compute_options["save_aug"] },
  aug_thinning { compute_options["aug_thinning"] },
  swap_leap { compute_options["swap_leap"] } ,
  log_aug_prob { zeros(dat.n_assessors) } {
    if(dat.any_missing){
      set_up_missing(dat.rankings, missing_indicator);
      initialize_missing_ranks(dat.rankings, missing_indicator);
    }
    if(save_aug){
      unsigned int nmc{ compute_options["nmc"] };
      augmented_data.set_size(dat.n_items, dat.n_assessors,
                              std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
      augmented_data.slice(0) = dat.rankings;
    }}

void Augmentation::augment_pairwise(
    const unsigned int t,
    Data& dat,
    const Parameters& pars,
    const Clustering& clus,
    const std::unique_ptr<Distance>& distfun
){
  if(!augpair) return;
  for(size_t i = 0; i < dat.n_assessors; ++i) {
    vec proposal;
    int g_diff{};
    if(pars.error_model == "none"){
      proposal = propose_pairwise_augmentation(
        dat.rankings.col(i), dat.items_above[i], dat.items_below[i]);
    } else if(pars.error_model == "bernoulli"){
      proposal = propose_swap(dat.rankings.col(i), dat.items_above[i],
                              dat.items_below[i], g_diff, swap_leap);
    } else {
      Rcpp::stop("error_model must be 'none' or 'bernoulli'");
    }

    double u = std::log(R::runif(0, 1));
    int cluster = clus.current_cluster_assignment(i);

    const vec& rankings = dat.rankings.col(i);
    double newdist = distfun->d(proposal, pars.rho_old.col(cluster));
    double olddist = distfun->d(rankings, pars.rho_old.col(cluster));
    double ratio = -pars.alpha_old(cluster) / dat.n_items * (newdist - olddist);

    if(pars.error_model != "none") {
      ratio += g_diff * std::log(pars.theta(t) / (1 - pars.theta(t)));
    }

    if(ratio > u) dat.rankings.col(i) = proposal;
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
      distfun, log_aug_prob(i)
    );
  }
}


