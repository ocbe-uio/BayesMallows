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

  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(pars.metric);
  auto rho_proposal = choose_rank_proposal(
    pars.rho_proposal_option, pars.leap_size);
  std::unique_ptr<ProposalDistribution> aug_prop;
  if(pars.error_model == "none"){
    aug_prop = std::make_unique<LeapAndShift>(1);
  } else if(pars.error_model == "bernoulli"){
    aug_prop = std::make_unique<Swap>(aug.swap_leap);
  } else {
    Rcpp::stop("error_model must be 'none' or 'bernoulli'");
  }

  clus.update_dist_mat(dat, pars, distfun);
  int aug_index = 0, cluster_assignment_index = 0;

  for(pars.t = 1; pars.t < pars.nmc; pars.t++){
    if (pars.t % 1000 == 0) {
      Rcpp::checkUserInterrupt();
      if(verbose){
        Rcpp::Rcout << "First " << pars.t <<
          " iterations of Metropolis-Hastings algorithm completed." << std::endl;
      }
    }

    pars.update_shape(dat, pris);
    pars.update_rho(dat, clus.current_cluster_assignment,
                    distfun, rho_proposal);
    pars.update_alpha(dat, distfun, pfun, pris,
                      clus.current_cluster_assignment);


  if(clus.clustering){
    clus.update_cluster_probs(pars, pris);
    clus.update_cluster_labels(pars.t, dat, pars, pfun);

    if(pars.t % clus.clus_thinning == 0){
      ++cluster_assignment_index;
      clus.cluster_assignment.col(cluster_assignment_index) = clus.current_cluster_assignment;
      clus.cluster_probs.col(cluster_assignment_index) = clus.current_cluster_probs;
    }
  }

  clus.update_wcd(pars.t);
  aug.update_missing_ranks(dat, clus, pars, distfun);
  aug.augment_pairwise(dat, pars, clus, distfun, aug_prop);

  if(aug.save_aug & (pars.t % aug.aug_thinning == 0)){
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
    Rcpp::Named("augmented_data") = aug.augmented_data
  );

}
