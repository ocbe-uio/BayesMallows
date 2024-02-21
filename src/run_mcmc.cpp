#include <RcppArmadillo.h>
#include "classes.h"
#include "rank_proposal.h"

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
  std::string metric = model_options["metric"];

  auto pfun = choose_partition_function(
    dat.n_items, metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(metric);
  auto rho_proposal = choose_rank_proposal(
    pars.rho_proposal_option, pars.leap_size);

  clus.update_dist_mat(dat, pars, distfun);
  int alpha_index = 0, rho_index = 0, aug_index = 0,
    cluster_assignment_index = 0;

  for(size_t t{1}; t < pars.nmc; ++t){
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t <<
          " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    pars.update_shape(t, dat, pris);
    pars.update_rho(t, rho_index, dat, clus.current_cluster_assignment,
                    distfun, rho_proposal);

    if(t % pars.alpha_jump == 0) {
      ++alpha_index;
      pars.update_alpha(alpha_index, dat, distfun, pfun, pris,
                        clus.current_cluster_assignment);
      pars.alpha_old = pars.alpha.col(alpha_index);
    }

  if(clus.clustering){
    clus.update_cluster_probs(pars, pris);
    clus.update_cluster_labels(t, dat, pars, pfun);

    if(t % clus.clus_thinning == 0){
      ++cluster_assignment_index;
      clus.cluster_assignment.col(cluster_assignment_index) = clus.current_cluster_assignment;
      clus.cluster_probs.col(cluster_assignment_index) = clus.current_cluster_probs;
    }
  }

  clus.update_wcd(t);
  aug.update_missing_ranks(dat, clus, pars, distfun);
  aug.augment_pairwise(t, dat, pars, clus, distfun);

  if(aug.save_aug & (t % aug.aug_thinning == 0)){
    ++aug_index;
    aug.augmented_data.slice(aug_index) = dat.rankings;
  }

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
    Rcpp::Named("any_missing") = dat.any_missing,
    Rcpp::Named("augpair") = aug.augpair,
    Rcpp::Named("n_assessors") = dat.n_assessors,
    Rcpp::Named("observation_frequency") = dat.observation_frequency
  );

}
