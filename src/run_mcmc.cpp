#include <RcppArmadillo.h>
#include "misc.h"
#include "mixtures.h"
#include "distances.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "parameterupdates.h"

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

  // The number of items ranked
  mat rankings = data["rankings"];
  rankings = rankings.t();
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  Rcpp::Nullable<vec> cardinalities = logz_list["cardinalities"];
  Rcpp::Nullable<vec> logz_estimate = logz_list["logz_estimate"];

  Rcpp::List constraints = data["constraints"];
  bool augpair = (constraints.length() > 0);
  bool any_missing = !is_finite(rankings);

  umat missing_indicator{};
  cube augmented_data{};
  bool save_aug = compute_options["save_aug"];
  int aug_thinning = compute_options["aug_thinning"];
  int rho_thinning = compute_options["rho_thinning"];
  int nmc = compute_options["nmc"];
  int n_clusters = model["n_clusters"];

  if(any_missing){
    set_up_missing(rankings, missing_indicator);
    initialize_missing_ranks(rankings, missing_indicator);
  }

  if(save_aug){
    augmented_data.set_size(n_items, n_assessors, std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
    augmented_data.slice(0) = rankings;
  }

  // Declare the cube to hold the latent ranks

  cube rho(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
  Rcpp::Nullable<mat> rho_init = initial_values["rho_init"];
  rho.slice(0) = initialize_rho(n_items, n_clusters, rho_init);
  mat rho_old = rho(span::all, span::all, span(0));

  // Declare the vector to hold the scaling parameter alpha
  int alpha_jump = compute_options["alpha_jump"];
  mat alpha(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
  double alpha_init = initial_values["alpha_init"];
  alpha.col(0).fill(alpha_init);

  // Clustering
  bool clustering = n_clusters > 1;
  int clus_thinning = compute_options["clus_thinning"];
  int n_cluster_assignments = n_clusters > 1 ? std::ceil(static_cast<double>(nmc * 1.0 / clus_thinning)) : 1;
  mat cluster_probs(n_clusters, n_cluster_assignments);
  cluster_probs.col(0).fill(1.0 / n_clusters);
  vec current_cluster_probs = cluster_probs.col(0);
  umat cluster_assignment(n_assessors, n_cluster_assignments);
  cluster_assignment.col(0) = randi<uvec>(n_assessors, distr_param(0, n_clusters - 1));
  uvec current_cluster_assignment = cluster_assignment.col(0);

  // Matrix with precomputed distances d(R_j, \rho_j), used to avoid looping during cluster assignment
  mat dist_mat(n_assessors, n_clusters);
  std::string metric = model["metric"];
  vec observation_frequency = data["observation_frequency"];
  update_dist_mat(dist_mat, rankings, rho_old, metric, observation_frequency);
  bool include_wcd = compute_options["include_wcd"];

  mat within_cluster_distance(n_clusters, include_wcd ? nmc : 1);
  within_cluster_distance.col(0) = update_wcd(current_cluster_assignment, dist_mat);


  int kappa_1 = priors["kappa_1"];
  int kappa_2 = priors["kappa_2"];

  // Declare vector with Bernoulli parameter for the case of intransitive preferences
  vec theta, shape_1, shape_2;
  std::string error_model = model["error_model"];
  if(error_model == "bernoulli"){
    theta = zeros<vec>(nmc);
    shape_1 = zeros<vec>(nmc);
    shape_2 = zeros<vec>(nmc);
    shape_1(0) = kappa_1;
    shape_2(0) = kappa_2;
  } else {
    theta.reset();
    shape_1.reset();
    shape_2.reset();
  }

  // Other variables used
  int alpha_index = 0, rho_index = 0, aug_index = 0, cluster_assignment_index = 0;
  vec alpha_old = alpha.col(0);

  uvec element_indices = regspace<uvec>(0, rankings.n_rows - 1);

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(int t = 1; t < nmc; ++t){
    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    if(error_model == "bernoulli"){
      update_shape_bernoulli(shape_1(t), shape_2(t), kappa_1, kappa_2,
                             rankings, constraints);

      // Update the theta parameter for the error model, which is independent of cluster
      theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
    }

    for(int i = 0; i < n_clusters; ++i){
      int leap_size = compute_options["leap_size"];
      update_rho(rho, rho_old, rho_index, i,
                 rho_thinning, alpha_old(i), leap_size,
                 clustering ? rankings.submat(element_indices, find(current_cluster_assignment == i)) : rankings,
                 metric, t, element_indices, observation_frequency);
    }

    if(t % alpha_jump == 0) {
      ++alpha_index;
      for(int i = 0; i < n_clusters; ++i){
        double lambda = priors["lambda"];
        double alpha_max = priors["alpha_max"];
        double alpha_prop_sd = compute_options["alpha_prop_sd"];

        alpha(i, alpha_index) = update_alpha(alpha_old(i),
              clustering ? rankings.submat(element_indices, find(current_cluster_assignment == i)) : rankings,
              clustering ? observation_frequency(find(current_cluster_assignment == i)) : observation_frequency,
              rho_old.col(i), alpha_prop_sd, metric, lambda, cardinalities, logz_estimate);

      }
      // Update alpha_old
      alpha_old = alpha.col(alpha_index);
    }

  if(clustering){
    bool save_ind_clus = compute_options["save_ind_clus"];
    int psi = priors["psi"];
    current_cluster_probs = update_cluster_probs(current_cluster_assignment, n_clusters, psi);

    current_cluster_assignment = update_cluster_labels(dist_mat, current_cluster_probs,
                                                       alpha_old, n_items, t, metric, cardinalities,
                                                       logz_estimate, save_ind_clus);

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
  if(any_missing){
    update_missing_ranks(rankings, current_cluster_assignment, missing_indicator,
                         alpha_old, rho_old, metric);
  }

  // Perform data augmentation of pairwise comparisons, if needed
  if(augpair){
    int swap_leap = compute_options["swap_leap"];
    augment_pairwise(rankings, current_cluster_assignment, alpha_old, 0.1, rho_old,
                     metric, constraints, error_model, swap_leap);
  }

  // Save augmented data if the user wants this. Uses the same index as rho.
  if(save_aug & (t % aug_thinning == 0)){
    ++aug_index;
    augmented_data.slice(aug_index) = rankings;
  }

  if(clustering | include_wcd){
    update_dist_mat(dist_mat, rankings, rho_old, metric, observation_frequency);
    }
  }

  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("shape1") = shape_1,
    Rcpp::Named("shape2") = shape_2,
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = cluster_probs,
    Rcpp::Named("within_cluster_distance") = within_cluster_distance,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("observation_frequency") = observation_frequency
  );
}
