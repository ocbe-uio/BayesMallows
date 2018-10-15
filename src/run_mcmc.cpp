#include <cmath>
#include "RcppArmadillo.h"
#include "clustering.h"
#include "distfuns.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"
#include "parameterupdates.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Worker function for computing the posterior distribution.
//'
//' @param rankings A set of complete rankings, with one sample per column.
//' With n_assessors samples and n_items items, rankings is n_items x n_assessors.
//' @param nmc Number of Monte Carlo samples.
//' @param constraints List of lists of lists, returned from `generate_constraints`.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function. Defaults to
//' \code{R_NilValue}.
//' @param logz_estimate Estimate of the log partition function.
//' @param metric The distance metric to use. One of \code{"spearman"},
//' \code{"footrule"}, \code{"kendall"}, \code{"cayley"}, or
//' \code{"hamming"}.
//' @param n_clusters Number of clusters. Defaults to 1.
//' @param include_wcd Boolean defining whether or
//' not to store the within-cluster distance.
//' @param leap_size Leap-and-shift step size.
//' @param alpha_prop_sd Standard deviation of proposal distribution for alpha.
//' @param alpha_init Initial value of alpha.
//' @param alpha_jump How many times should we sample \code{rho} between
//' each time we sample \code{alpha}. Setting \code{alpha_jump} to a high
//' number can significantly speed up computation time, since we then do not
//' have to do expensive computation of the partition function.
//' @param lambda Parameter of the prior distribution.
//' @param psi Hyperparameter for the Dirichlet prior distribution used in clustering.
//' @param rho_thinning Thinning parameter. Keep only every \code{rho_thinning} rank
//' sample from the posterior distribution.
//' @param aug_thinning Integer specifying the thinning for data augmentation.
//' @param cluster_assignment_thinning Integer specifying the thinning for saving cluster assignments.
//' @param save_augmented_data Whether or not to save the augmented data every
//' \code{aug_thinning}th iteration.
//' @param verbose Logical specifying whether to print out the progress of the
//' Metropolis-Hastings algorithm. If \code{TRUE}, a notification is printed every
//' 1000th iteration.
//' @keywords internal
//'
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat rankings, int nmc,
                    Rcpp::List constraints,
                    Rcpp::Nullable<arma::vec> cardinalities,
                    Rcpp::Nullable<arma::vec> logz_estimate,
                    Rcpp::Nullable<arma::vec> rho_init,
                    std::string metric = "footrule",
                    int n_clusters = 1,
                    bool include_wcd = false,
                    int leap_size = 1,
                    double alpha_prop_sd = 0.5,
                    double alpha_init = 5,
                    int alpha_jump = 1,
                    double lambda = 0.1,
                    int psi = 10,
                    int rho_thinning = 1,
                    int aug_thinning = 1,
                    int cluster_assignment_thinning = 1,
                    bool save_augmented_data = false,
                    bool verbose = false
                      ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  // Number of alpha values to store, per cluster.
  int n_alpha = ceil(nmc * 1.0 / alpha_jump);

  // Number of rho values to store, per cluster and item
  int n_rho = ceil(nmc * 1.0 / rho_thinning);

  // Number of augmented data sets to store
  int n_aug = ceil(nmc * 1.0 / aug_thinning);

  // Number of cluster assignments to store
  int n_cluster_assignments = ceil(nmc * 1.0 / cluster_assignment_thinning);

  // Check if we want to do clustering
  bool clustering = n_clusters > 1;

  // Check if we have pairwise preferences
  bool augpair;

  if(constraints.length() > 0){
    augpair = true;
  } else {
    augpair = false;
  }

  // Boolean which indicates if ANY assessor has missing ranks
  bool any_missing = !arma::is_finite(rankings);

  arma::mat missing_indicator;
  arma::vec assessor_missing;

  if(any_missing){
    missing_indicator = arma::zeros<arma::mat>(n_items, n_assessors);

    // Number of missing items per assessor
    assessor_missing = arma::zeros<arma::vec>(n_assessors);

    // Fill the two above defined missingness indicators
    define_missingness(missing_indicator, assessor_missing, rankings, n_items, n_assessors);

  } else {
    missing_indicator.reset();
    assessor_missing.reset();
  }

  // Declare the matrix to hold the latent ranks
  // Note: Armadillo matrices are stored in column-major ordering. Hence,
  // we put the items along the column, since they are going to be accessed at the
  // same time for a given Monte Carlo sample.
  arma::cube rho(n_items, n_clusters, n_rho);

  // Declare the vector to hold the scaling parameter alpha
  arma::mat alpha(n_clusters, n_alpha);

  // Set the initial alpha value
  alpha.col(0).fill(alpha_init);

  // Initialize latent ranks as provided by rho_init, or randomly:
  for(int i = 0; i < n_clusters; ++i){
    if(rho_init.isNotNull()){
      rho.slice(0).col(i) = Rcpp::as<arma::vec>(rho_init);
    } else {
      rho.slice(0).col(i) = arma::shuffle(arma::regspace<arma::vec>(1, 1, n_items));
    }
  }

  // Fill in missing ranks, if needed
  if(any_missing){
    initialize_missing_ranks(rankings, missing_indicator, assessor_missing,
                             n_items, n_assessors);
  }


  // If the user wants to save augmented data, we need a cube
  arma::cube augmented_data;
  if(save_augmented_data){
    augmented_data.set_size(n_items, n_assessors, n_aug);
    augmented_data.slice(0) = rankings;
  }

  // Cluster probabilities
  arma::mat cluster_probs;

  // Declare the cluster indicator z
  arma::umat cluster_assignment;
  arma::uvec current_cluster_assignment;

  // Within cluster distance
  arma::mat within_cluster_distance;

  // Submatrix of rankings to be updated in each clustering step
  arma::mat clus_mat = rankings;

  // Matrix with precomputed distances d(R_j, \rho_j), used to avoid looping during cluster assignment
  arma::mat dist_mat;

  if(clustering | include_wcd){

    cluster_assignment.set_size(n_assessors, n_cluster_assignments);
    within_cluster_distance.set_size(n_clusters, nmc);

    if(clustering){
      cluster_probs.set_size(n_clusters, nmc);
      cluster_probs.col(0).fill(1.0/n_clusters);

      // Initialize clusters randomly
      cluster_assignment.col(0) = arma::randi<arma::uvec>(n_assessors, arma::distr_param(0, n_clusters - 1));

    } else {
      // Set all clusters once and for all
      cluster_assignment = arma::zeros<arma::umat>(n_assessors, 1);
    }

    // Initialize the distance matrix. Can be done here since \rho and R already are intialized
    dist_mat.set_size(n_assessors, n_clusters);

    for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){
      update_distance_matrix(dist_mat, rankings, rho.slice(0).col(cluster_index),
                             n_assessors, cluster_index, metric);
    }

    current_cluster_assignment = cluster_assignment.col(0);

    update_wcd(within_cluster_distance, current_cluster_assignment,
               dist_mat, n_clusters, 0);

  }


  // Declare indicator vectors to hold acceptance or not
  arma::vec alpha_acceptance = arma::ones(n_clusters);
  arma::vec rho_acceptance = arma::ones(n_clusters);

  arma::vec aug_acceptance;
  if(any_missing | augpair){
    aug_acceptance = arma::ones<arma::vec>(n_assessors);
  } else {
    aug_acceptance.reset();
  }

  // Other variables used
  int alpha_index = 0, rho_index = 0, aug_index = 0, cluster_assignment_index = 0;
  arma::vec alpha_old = alpha.col(0);
  arma::mat rho_old = rho(arma::span::all, arma::span::all, arma::span(0));
  bool rho_accepted = false;
  bool augmentation_accepted = false;

  arma::uvec element_indices = arma::regspace<arma::uvec>(0, rankings.n_rows - 1);

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(int t = 1; t < nmc; ++t){

    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << t
        << " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }

    }

    if(clustering){
      update_cluster_probs(cluster_probs, current_cluster_assignment, n_clusters, psi, t);
    }


    for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){

      // Find the members of this cluster
      arma::uvec matches = arma::find(current_cluster_assignment == cluster_index);

      // Matrix of ranks for this cluster
      if(clustering){
        clus_mat = rankings.submat(element_indices, matches);
      } else if (any_missing | augpair){
        // When augmenting data, we need to update in every step, even without clustering
        clus_mat = rankings;
      }

      // Call the void function which updates rho by reference
      update_rho(rho, rho_acceptance, rho_old, rho_index, cluster_index,
                 rho_thinning, alpha_old(cluster_index), leap_size, clus_mat, metric, n_items, t,
                 element_indices, rho_accepted);

      if((rho_accepted | augmentation_accepted) & (clustering | include_wcd)){
        // Note: Must use rho_old rather than rho, because when rho_thinning > 1,
        // rho does not necessarily have the last accepted value
        update_distance_matrix(dist_mat, rankings, rho_old.col(cluster_index),
                               n_assessors, cluster_index, metric);
      }

      if(t % alpha_jump == 0) {

        // Increment alpha_index only once, and not n_cluster times!
        if(cluster_index == 0) ++alpha_index;

        // Call the void function which updates alpha by reference
        update_alpha(alpha, alpha_acceptance, alpha_old, clus_mat, alpha_index,
                     cluster_index, rho_old, alpha_prop_sd, metric, lambda, n_items,
                     cardinalities, logz_estimate);
      }

    }


  if(clustering){
    // Update the cluster labels, per assessor
    update_cluster_labels(cluster_assignment, current_cluster_assignment, dist_mat, rho_old, rankings, cluster_probs,
                          alpha_old, n_items, n_assessors, n_clusters, cluster_assignment_thinning, cluster_assignment_index,
                          t, metric, cardinalities, logz_estimate);
  }

  if(include_wcd){
    // Update within_cluster_distance
    update_wcd(within_cluster_distance, current_cluster_assignment,
               dist_mat, n_clusters, t);
  }

  // Perform data augmentation of missing ranks, if needed
  if(any_missing){
    update_missing_ranks(rankings, current_cluster_assignment, aug_acceptance, missing_indicator,
                         assessor_missing, n_items, n_assessors, alpha_old, rho_old,
                         metric, t, clustering, augmentation_accepted);
  }

    // Perform data augmentation of pairwise comparisons, if needed
  if(augpair){
    augment_pairwise(rankings, current_cluster_assignment, alpha_old, rho_old,
                     metric, constraints, n_assessors, n_items, t,
                     aug_acceptance, clustering, augmentation_accepted);

  }

  // Save augmented data if the user wants this. Uses the same index as rho.
  if(save_augmented_data & (t % aug_thinning == 0)){
    ++aug_index;
    augmented_data.slice(aug_index) = rankings;
  }
  }


  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance / nmc,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance / nmc,
    Rcpp::Named("cluster_assignment") = cluster_assignment + 1,
    Rcpp::Named("cluster_probs") = cluster_probs,
    Rcpp::Named("within_cluster_distance") = within_cluster_distance,
    Rcpp::Named("augmented_data") = augmented_data,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("aug_acceptance") = aug_acceptance / nmc,
    Rcpp::Named("n_assessors") = n_assessors
  );


}



