#include <math.h>
#include "RcppArmadillo.h"
#include "parameterupdates.h"
#include "misc.h"
#include "missing_data.h"
#include "pairwise_comparisons.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

//' Worker function for computing the posterior distribtuion.
//'
//' @param rankings A set of complete rankings, with one sample per column.
//' With n_assessors samples and n_items items, rankings is n_items x n_assessors.
//' @param nmc Number of Monte Carlo samples.
//' @param preferences Matrix of preferences preferences, 3 x n_items.
//' @param constrained Matrix of constrained elements, 2 rows.
//' @param cardinalities Used when metric equals \code{"footrule"} or
//' \code{"spearman"} for computing the partition function. Defaults to
//' \code{R_NilValue}.
//' @param is_fit Importance sampling fit.
//' @param metric The distance metric to use. One of \code{"spearman"},
//' \code{"footrule"}, \code{"kendall"}, \code{"cayley"}, or
//' \code{"hamming"}.
//' @param n_clusters Number of clusters. Defaults to 1.
//' @param leap_size Leap-and-shift step size.
//' @param sd_alpha Standard deviation of proposal distribution for alpha.
//' @param alpha_init Initial value of alpha.
//' @param alpha_jump How many times should we sample \code{rho} between
//' each time we sample \code{alpha}. Setting \code{alpha_jump} to a high
//' number can significantly speed up computation time, since we then do not
//' have to do expensive computation of the partition function.
//' @param lambda Parameter of the prior distribution.
//' @param thinning Thinning parameter. Keep only every \code{thinning} rank
//' sample from the posterior distribution.
//' @param aug_diag_thinning The interval in which we save
//' augmentation diagnostics.
//'
// [[Rcpp::export]]
Rcpp::List run_mcmc(arma::mat rankings, int nmc,
                    Rcpp::Nullable<arma::mat> preferences,
                    Rcpp::Nullable<arma::mat> constrained,
                    Rcpp::Nullable<arma::vec> cardinalities,
                    Rcpp::Nullable<arma::vec> is_fit,
                    std::string metric = "footrule",
                    int n_clusters = 1,
                    int leap_size = 1,
                    double sd_alpha = 0.5,
                    double alpha_init = 5,
                    int alpha_jump = 1,
                    double lambda = 0.1,
                    int thinning = 1,
                    int aug_diag_thinning = 100
                      ){

  // The number of items ranked
  int n_items = rankings.n_rows;

  // The number of assessors
  int n_assessors = rankings.n_cols;

  // Number of alpha values to store, per cluster.
  int n_alpha = ceil(nmc * 1.0 / alpha_jump);

  // Number of rho values to store, per cluster and item
  int n_rho = ceil(nmc * 1.0 / thinning);

  // Number of augmentation diagnostics to store
  int n_aug_diag = ceil(nmc * 1.0 / aug_diag_thinning);

  // Check if we want to do clustering
  bool clustering = n_clusters > 1;

  // Check if we have pairwise preferences
  bool augpair;
  arma::mat pairwise_preferences, constrained_elements;


  if(preferences.isNotNull() & constrained.isNotNull()){
    augpair = true;
    pairwise_preferences = Rcpp::as<arma::mat>(preferences);
    constrained_elements = Rcpp::as<arma::mat>(constrained);
  } else {
    augpair = false;
  }

  // Declare indicator matrix of missing ranks, and fill it with zeros
  arma::mat missing_indicator = arma::zeros<arma::mat>(n_items, n_assessors);

  // Number of missing items per assessor
  arma::vec assessor_missing = arma::zeros<arma::vec>(n_assessors);

  // Fill the two above defined missingness indicators
  define_missingness(missing_indicator, assessor_missing, rankings, n_items, n_assessors);

  // Boolean which indicates if ANY assessor has missing ranks
  bool any_missing = arma::any(assessor_missing);

  // Declare the matrix to hold the latent ranks
  // Note: Armadillo matrices are stored in column-major ordering. Hence,
  // we put the items along the column, since they are going to be accessed at the
  // same time for a given Monte Carlo sample.
  arma::cube rho(n_items, n_clusters, n_rho);

  // Declare the vector to hold the scaling parameter alpha
  arma::mat alpha(n_clusters, n_alpha);

  // Set the initial alpha value
  alpha.col(0).fill(alpha_init);

  // Initialize latent ranks randomly
  for(int i = 0; i < n_clusters; ++i){
    rho.slice(0).col(i) = arma::shuffle(arma::regspace<arma::vec>(1, 1, n_items));
  }

  // Hyperparameter for Dirichlet prior used in clustering
  // Consider having this as an optional user argument
  int psi = 10;

  // Cluster probabilities
  arma::mat cluster_probs;
  // Declare the cluster indicator z
  arma::umat cluster_indicator;
  // Submatrix of rankings to be updated in each clustering step
  arma::mat clus_mat;

  if(clustering){
    cluster_probs.set_size(n_clusters, nmc);
    cluster_probs.col(0).fill(1.0/n_clusters);
    // Initialize clusters randomly
    cluster_indicator.set_size(n_assessors, nmc);
    cluster_indicator.col(0) = arma::randi<arma::uvec>(n_assessors, arma::distr_param(0, n_clusters - 1));

  } else {
    // Empty the matrix
    cluster_probs.reset();
    cluster_indicator.reset();
    clus_mat = rankings; // assign once and for all
  }




  // Fill in missing ranks, if needed
  if(any_missing){
    initialize_missing_ranks(rankings, missing_indicator, assessor_missing,
                             n_items, n_assessors);
  }

  // Declare indicator vectors to hold acceptance or not
  arma::vec alpha_acceptance = arma::ones(n_clusters);
  arma::vec rho_acceptance = arma::ones(n_clusters);
  arma::mat aug_acceptance = arma::zeros<arma::mat>(n_assessors, n_aug_diag);

  // Other variables used
  int alpha_index = 0, rho_index = 0, aug_diag_index = 0;
  arma::vec alpha_old = alpha.col(0);
  arma::mat rho_old = rho(arma::span::all, arma::span::all, arma::span(0));

  arma::uvec element_indices = arma::regspace<arma::uvec>(0, rankings.n_rows - 1);

  // This is the Metropolis-Hastings loop

  // Starting at t = 1, meaning that alpha and rho must be initialized at index 0,
  // and this has been done above
  for(int t = 1; t < nmc; ++t){

    // Check if the user has tried to interrupt.
    if (t % 1000 == 0) Rcpp::checkUserInterrupt();

    if(clustering){
      update_cluster_probs(cluster_probs, cluster_indicator, n_clusters, psi, t);
    }


    for(int cluster_index = 0; cluster_index < n_clusters; ++cluster_index){

      // Matrix of ranks for this cluster
      if(clustering){
        clus_mat = rankings.submat(element_indices,
                                   arma::find(cluster_indicator.col(t - 1) == cluster_index));
      }


      // Call the void function which updates rho by reference
      update_rho(rho, rho_acceptance, rho_old, rho_index, cluster_index,
                 thinning, alpha_old(cluster_index), leap_size, clus_mat, metric, n_items, t,
                 element_indices);

      if(t % alpha_jump == 0) {

        // Increment alpha_index only once, and not n_cluster times!
        if(cluster_index == 0) ++alpha_index;


        // Call the void function which updates alpha by reference
        update_alpha(
          alpha, alpha_acceptance, alpha_old,
          clus_mat,
          alpha_index, cluster_index, rho_old,
          sd_alpha, metric, lambda, n_items,
          cardinalities, is_fit);
      }

    }

  // Update the cluster labels, per assessor
  if(clustering){
    update_cluster_labels(cluster_indicator, rho_old, rankings, cluster_probs,
                          alpha_old, n_items, n_assessors, n_clusters,
                          t, metric, cardinalities, is_fit);
  }



  // Perform data augmentation of missing ranks, if needed
  if(any_missing){
    update_missing_ranks(rankings, cluster_indicator, aug_acceptance, missing_indicator,
                         assessor_missing, n_items, n_assessors, alpha_old, rho_old,
                         metric, t, aug_diag_index, aug_diag_thinning, clustering);
  }


    // Perform data augmentation of pairwise comparisons, if needed
  if(augpair){
    augment_pairwise(rankings, cluster_indicator, alpha_old, rho_old,
                     metric, pairwise_preferences, constrained_elements,
                     n_assessors, n_items, t, aug_acceptance, aug_diag_index,
                     aug_diag_thinning);
  }


  }


  // Return everything that might be of interest
  return Rcpp::List::create(
    Rcpp::Named("rho") = rho,
    Rcpp::Named("rho_acceptance") = rho_acceptance/nmc,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("alpha_acceptance") = alpha_acceptance/nmc,
    Rcpp::Named("cluster_indicator") = cluster_indicator + 1,
    Rcpp::Named("cluster_probs") = cluster_probs,
    Rcpp::Named("any_missing") = any_missing,
    Rcpp::Named("augpair") = augpair,
    Rcpp::Named("aug_acceptance") = aug_acceptance,
    Rcpp::Named("metric") = metric,
    Rcpp::Named("lambda") = lambda,
    Rcpp::Named("nmc") = nmc,
    Rcpp::Named("n_items") = n_items,
    Rcpp::Named("n_assessors") = n_assessors,
    Rcpp::Named("n_clusters") = n_clusters,
    Rcpp::Named("alpha_jump") = alpha_jump,
    Rcpp::Named("thinning") = thinning,
    Rcpp::Named("leap_size") = leap_size,
    Rcpp::Named("sd_alpha") = sd_alpha,
    Rcpp::Named("aug_diag_thinning") = aug_diag_thinning
  );

// return(Rcpp::List::create(Rcpp::Named("test") = 0));
}



