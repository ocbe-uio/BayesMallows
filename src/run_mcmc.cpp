#include <RcppArmadillo.h>
#include "classes.h"
#include "rank_proposal.h"
#include "progress_reporter.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc(
    Rcpp::List data,
    Rcpp::List model_options,
    Rcpp::List compute_options,
    Rcpp::List priors,
    Rcpp::List initial_values,
    Rcpp::Nullable<arma::mat> pfun_values,
    Rcpp::Nullable<arma::mat> pfun_estimate,
    bool verbose = false){
  Data dat{data};
  Priors pris{priors};
  Parameters pars{model_options, compute_options, initial_values, dat.n_items};
  Clustering clus{pars, compute_options, dat.n_assessors};
  Augmentation aug{dat, compute_options};
  ProgressReporter rep{verbose};

  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(pars.metric);
  auto rho_proposal = choose_rho_proposal(
    pars.rho_proposal_option, pars.leap_size);
  auto partial_aug_prop = choose_partial_proposal(aug.aug_method, aug.pseudo_aug_metric);
  auto pairwise_aug_prop = choose_pairwise_proposal(pars.error_model, aug.swap_leap);

  clus.update_dist_mat(dat, pars, distfun);

  for(pars.t = 1; pars.t < pars.nmc; pars.t++){
    rep.report(pars.t);
    pars.update_shape(dat, pris);
    pars.update_rho(dat, clus.current_cluster_assignment,
                    distfun, rho_proposal);
    pars.update_alpha(dat, distfun, pfun, pris,
                      clus.current_cluster_assignment);
    clus.update_cluster_probs(pars, pris);
    clus.update_cluster_labels(pars.t, dat, pars, pfun);
    clus.save_cluster_parameters(pars.t);
    clus.update_wcd(pars.t);
    aug.update_missing_ranks(dat, clus, pars, distfun, partial_aug_prop);
    aug.augment_pairwise(dat, pars, clus, distfun, pairwise_aug_prop);
    aug.save_augmented_data(dat, pars);
    clus.update_dist_mat(dat, pars, distfun);
  }

  return Rcpp::List::create(
    Rcpp::Named("rho") = pars.rho,
    Rcpp::Named("alpha") = pars.alpha,
    Rcpp::Named("theta") = pars.theta,
    Rcpp::Named("shape1") = pars.shape_1,
    Rcpp::Named("shape2") = pars.shape_2,
    Rcpp::Named("cluster_assignment") = clus.cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = clus.cluster_probs,
    Rcpp::Named("within_cluster_distance") = clus.within_cluster_distance,
    Rcpp::Named("augmented_data") = aug.augmented_data,
    Rcpp::Named("alpha_acceptance") = pars.alpha_acceptance /
      (pars.nmc - pars.burnin) * pars.alpha_jump,
    Rcpp::Named("rho_acceptance") = pars.rho_acceptance / (pars.nmc - pars.burnin)
  );

}
