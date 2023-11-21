#include <RcppArmadillo.h>
#include "misc.h"
#include "mixtures.h"
#include "distances.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "parameterupdates.h"
#include "parameters.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List run_mcmc(Rcpp::List data,
                    Rcpp::List model,
                    Rcpp::List compute_options,
                    Rcpp::List priors,
                    Rcpp::List initial_values,
                    Rcpp::List logz_list,
                    bool verbose = false
                      ){


  Data dat{data, compute_options};
  Priors pris{priors};
  Parameters pars{model, compute_options, initial_values, dat.n_items};
  Clustering clus{pars, compute_options, dat.n_assessors};

  update_dist_mat(clus.dist_mat, dat.rankings, pars.rho_old, pars.metric, dat.observation_frequency);

  int alpha_index = 0, rho_index = 0, aug_index = 0, cluster_assignment_index = 0;

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(size_t t{1}; t < pars.nmc; ++t){
    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    pars.update_shape(t, dat, pris);

    for(int i = 0; i < pars.n_clusters; ++i){
      pars.update_rho(i, t, rho_index, dat);
    }

    if(t % pars.get_alpha_jump() == 0) {
      ++alpha_index;
      for(int i = 0; i < pars.n_clusters; ++i){
        pars.update_alpha(i, alpha_index, dat, logz_list, pris);
      }
      // Update alpha_old
      pars.alpha_old = pars.alpha.col(alpha_index);
    }

  if(clus.clustering){

    clus.update_cluster_probs(pars, pris);
    clus.update_cluster_labels(t, dat, pars, logz_list);

    if(t % clus.clus_thinning == 0){
      ++cluster_assignment_index;
      clus.cluster_assignment.col(cluster_assignment_index) = clus.current_cluster_assignment;
      clus.cluster_probs.col(cluster_assignment_index) = clus.current_cluster_probs;
    }
  }

  clus.update_wcd(t);

  // Perform data augmentation of missing ranks, if needed
  if(dat.any_missing){
    update_missing_ranks(dat.rankings, clus.current_cluster_assignment, dat.missing_indicator,
                         pars.alpha_old, pars.rho_old, pars.metric);
  }

  // Perform data augmentation of pairwise comparisons, if needed
  if(dat.augpair){
    int swap_leap = compute_options["swap_leap"];
    augment_pairwise(dat.rankings, clus.current_cluster_assignment, pars.alpha_old, 0.1, pars.rho_old,
                     pars.metric, dat.constraints, pars.get_error_model(), swap_leap);
  }

  // Save augmented data if the user wants this. Uses the same index as rho.
  if(dat.save_aug & (t % dat.aug_thinning == 0)){
    ++aug_index;
    dat.augmented_data.slice(aug_index) = dat.rankings;
  }

  if(clus.clustering | clus.include_wcd){
    update_dist_mat(clus.dist_mat, dat.rankings, pars.rho_old, pars.metric, dat.observation_frequency);
    }
  }

  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = pars.rho,
    Rcpp::Named("alpha") = pars.alpha,
    Rcpp::Named("theta") = pars.theta,
    Rcpp::Named("shape1") = pars.shape_1,
    Rcpp::Named("shape2") = pars.shape_2,
    Rcpp::Named("cluster_assignment") = clus.cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = clus.cluster_probs,
    Rcpp::Named("within_cluster_distance") = clus.within_cluster_distance,
    Rcpp::Named("augmented_data") = dat.augmented_data,
    Rcpp::Named("any_missing") = dat.any_missing,
    Rcpp::Named("augpair") = dat.augpair,
    Rcpp::Named("n_assessors") = dat.n_assessors,
    Rcpp::Named("observation_frequency") = dat.observation_frequency
  );
}
