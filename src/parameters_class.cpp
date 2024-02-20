#include "classes.h"
#include "proposal_functions.h"
#include "distributions.h"
using namespace arma;

Parameters::Parameters(
  const Rcpp::List& model_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values,
  const unsigned int n_items) :
  n_clusters { model_options["n_clusters"] },
  nmc { compute_options["nmc"] },
  error_model(model_options["error_model"]),
  alpha_jump { compute_options["alpha_jump"] },
  leap_size { compute_options["leap_size"] },
  rho_proposal_option( compute_options["rho_proposal"] ),
  element_indices { regspace<uvec>(0, n_items - 1) },
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  rho_thinning { compute_options["rho_thinning"] }
  {
    alpha.set_size(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
    double alpha_init = initial_values["alpha_init"];
    alpha.col(0).fill(alpha_init);
    alpha_old = alpha.col(0);

    rho.set_size(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
    Rcpp::Nullable<mat> rho_init = initial_values["rho_init"];
    if(rho_init.isNotNull()){
      rho.slice(0) = repmat(Rcpp::as<mat>(rho_init), 1, n_clusters);
    } else {
      for (size_t i{}; i < n_clusters; ++i) {
        ivec a = Rcpp::sample(n_items, n_items);
        rho.slice(0).col(i) = conv_to<vec>::from(a);
      }
    }
    rho_old = rho(span::all, span::all, span(0));

    if(error_model == "bernoulli"){
      theta = zeros<vec>(nmc);
      shape_1 = zeros<vec>(nmc);
      shape_2 = zeros<vec>(nmc);
    } else {
      theta.reset();
      shape_1.reset();
      shape_2.reset();
    }
  }

void Parameters::update_rho(
    int t,
    int& rho_index,
    const Data& dat,
    const uvec& current_cluster_assignment,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<ProposalDistribution>& prop
) {
  for(size_t i{}; i < n_clusters; ++i){
    const uvec cluster_indicator = find(current_cluster_assignment == i);
    const mat cluster_rankings = dat.rankings.submat(
      element_indices, cluster_indicator);
    const vec cluster_frequency =
      dat.observation_frequency.elem(cluster_indicator);
    vec rho_cluster = rho_old.col(i);
    rho_old.col(i) = make_new_rho(rho_cluster, cluster_rankings, alpha_old(i),
                leap_size, distfun, cluster_frequency);
    if(t % rho_thinning == 0){
      if(i == 0) ++rho_index;
      rho.slice(rho_index).col(i) = rho_old.col(i);
    }
  }
}

void Parameters::update_shape(
    int t, const Data& dat, const Priors& priors) {
  if(error_model != "bernoulli") return;
  int sum_1{};
  int sum_2{};
  for(size_t i = 0; i < dat.n_assessors; ++i){
    for(size_t j = 0; j < dat.n_items; ++j) {
      uvec items_above = dat.items_above[i][j];
      uvec items_below = dat.items_below[i][j];

      for(size_t k = 0; k < items_above.n_elem; ++k){
        int g = (as_scalar(dat.rankings.col(i).row(j)) < as_scalar(dat.rankings.col(i).row(items_above(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
      for(size_t k = 0; k < items_below.n_elem; ++k){
        int g = (as_scalar(dat.rankings.col(i).row(j)) > as_scalar(dat.rankings.col(i).row(items_below(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
    }
  }
  shape_1(t) = priors.kappa(0) + sum_1;
  shape_2(t) = priors.kappa(1) + sum_2;
  theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
}

void Parameters::update_alpha(
    int alpha_index,
    const Data& dat,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartitionFunction>& pfun,
    const Priors& priors,
    const uvec& current_cluster_assignment) {

  for(size_t i{}; i < n_clusters; ++i) {
    const uvec cluster_indicator = find(current_cluster_assignment == i);
    const mat cluster_rankings = dat.rankings.submat(
      element_indices, cluster_indicator);
    const vec cluster_frequency =
      dat.observation_frequency.elem(cluster_indicator);

    AlphaRatio test = make_new_alpha(
      alpha_old(i), rho_old.col(i),
      alpha_prop_sd, distfun, pfun, cluster_rankings,
      cluster_frequency, dat.n_items, priors);

    if(test.accept){
      alpha(i, alpha_index) = test.proposal;
    } else {
      alpha(i, alpha_index) = alpha_old(i);
    }
  }
}
