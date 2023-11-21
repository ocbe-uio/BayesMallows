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

  // Clustering
  bool clustering = pars.get_n_clusters() > 1;
  int clus_thinning = compute_options["clus_thinning"];
  int n_cluster_assignments = pars.get_n_clusters() > 1 ? std::ceil(static_cast<double>(pars.get_nmc() * 1.0 / clus_thinning)) : 1;
  mat cluster_probs(pars.get_n_clusters(), n_cluster_assignments);
  cluster_probs.col(0).fill(1.0 / pars.get_n_clusters());
  vec current_cluster_probs = cluster_probs.col(0);
  umat cluster_assignment(dat.n_assessors, n_cluster_assignments);
  cluster_assignment.col(0) = randi<uvec>(dat.n_assessors, distr_param(0, pars.get_n_clusters() - 1));
  uvec current_cluster_assignment = cluster_assignment.col(0);

  // Matrix with precomputed distances d(R_j, \rho_j), used to avoid looping during cluster assignment
  mat dist_mat(dat.n_assessors, pars.get_n_clusters());
  vec observation_frequency = data["observation_frequency"];
  update_dist_mat(dist_mat, dat.rankings, pars.rho_old, pars.get_metric(), observation_frequency);
  bool include_wcd = compute_options["include_wcd"];

  mat within_cluster_distance(pars.get_n_clusters(), include_wcd ? pars.get_nmc() : 1);
  within_cluster_distance.col(0) = update_wcd(current_cluster_assignment, dist_mat);


  // Other variables used
  int alpha_index = 0, rho_index = 0, aug_index = 0, cluster_assignment_index = 0;

  uvec element_indices = regspace<uvec>(0, dat.rankings.n_rows - 1);

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(size_t t{1}; t < pars.get_nmc(); ++t){
    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    if(pars.get_error_model() == "bernoulli") pars.update_shape(t, dat.rankings, dat.constraints, pris);

    for(int i = 0; i < pars.get_n_clusters(); ++i){
      pars.update_rho(i, t, rho_index, dat.rankings, observation_frequency);
    }

    if(t % pars.get_alpha_jump() == 0) {
      ++alpha_index;
      for(int i = 0; i < pars.get_n_clusters(); ++i){
        pars.update_alpha(i, alpha_index, dat.rankings, observation_frequency, logz_list, pris);
      }
      // Update alpha_old
      pars.alpha_old = pars.alpha.col(alpha_index);
    }

  if(clustering){
    bool save_ind_clus = compute_options["save_ind_clus"];
    int psi = priors["psi"];
    current_cluster_probs = update_cluster_probs(current_cluster_assignment, pars.get_n_clusters(), psi);

    current_cluster_assignment = update_cluster_labels(
      dist_mat, current_cluster_probs, pars.alpha_old, dat.n_items, t, pars.get_metric(), logz_list, save_ind_clus);

    if(t % clus_thinning == 0){
      ++cluster_assignment_index;
      cluster_assignment.col(cluster_assignment_index) = current_cluster_assignment;
      cluster_probs.col(cluster_assignment_index) = current_cluster_probs;
    }
  }

  if(include_wcd){
    // Update within_cluster_distance
    within_cluster_distance.col(t) = update_wcd(current_cluster_assignment, dist_mat);
  }

  // Perform data augmentation of missing ranks, if needed
  if(dat.any_missing){
    update_missing_ranks(dat.rankings, current_cluster_assignment, dat.missing_indicator,
                         pars.alpha_old, pars.rho_old, pars.get_metric());
  }

  // Perform data augmentation of pairwise comparisons, if needed
  if(dat.augpair){
    int swap_leap = compute_options["swap_leap"];
    augment_pairwise(dat.rankings, current_cluster_assignment, pars.alpha_old, 0.1, pars.rho_old,
                     pars.get_metric(), dat.constraints, pars.get_error_model(), swap_leap);
  }

  // Save augmented data if the user wants this. Uses the same index as rho.
  if(dat.save_aug & (t % dat.aug_thinning == 0)){
    ++aug_index;
    dat.augmented_data.slice(aug_index) = dat.rankings;
  }

  if(clustering | include_wcd){
    update_dist_mat(dist_mat, dat.rankings, pars.rho_old, pars.get_metric(), observation_frequency);
    }
  }

  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = pars.rho,
    Rcpp::Named("alpha") = pars.alpha,
    Rcpp::Named("theta") = pars.theta,
    Rcpp::Named("shape1") = pars.shape_1,
    Rcpp::Named("shape2") = pars.shape_2,
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = cluster_probs,
    Rcpp::Named("within_cluster_distance") = within_cluster_distance,
    Rcpp::Named("augmented_data") = dat.augmented_data,
    Rcpp::Named("any_missing") = dat.any_missing,
    Rcpp::Named("augpair") = dat.augpair,
    Rcpp::Named("n_assessors") = dat.n_assessors,
    Rcpp::Named("observation_frequency") = observation_frequency
  );
}
