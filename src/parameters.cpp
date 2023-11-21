#include "parameters.h"

using namespace arma;

parameters::parameters(
  const Rcpp::List& model,
  const Rcpp::List& compute_options,
  const Rcpp::List& priors,
  const Rcpp::List& initial_values,
  const int n_items) :
  metric { verify_metric(model["metric"]) },
  n_items { n_items },
  nmc { verify_positive(compute_options["nmc"]) }
  {
    int n_clusters = model["n_clusters"];
    int alpha_jump = compute_options["alpha_jump"];
    alpha.set_size(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
    double alpha_init = initial_values["alpha_init"];
    alpha.col(0).fill(alpha_init);
    alpha_old = alpha.col(0);
    lambda = Rcpp::as<double>(priors["lambda"]);
    alpha_prop_sd = Rcpp::as<double>(compute_options["alpha_prop_sd"]);

    rho_thinning = Rcpp::as<int>(compute_options["rho_thinning"]);
    rho.set_size(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
    Rcpp::Nullable<mat> rho_init = initial_values["rho_init"];
    rho.slice(0) = initialize_rho(n_items, n_clusters, rho_init);
    rho_old = rho(span::all, span::all, span(0));

    leap_size = Rcpp::as<int>(compute_options["leap_size"]);

    kappa_1 = Rcpp::as<int>(priors["kappa_1"]);
    kappa_2 = Rcpp::as<int>(priors["kappa_2"]);

    error_model = Rcpp::as<std::string>(model["error_model"]);
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
  }

void parameters::update_rho(int cluster_index, int t, int& rho_index,
                            const mat& rankings,
                            const vec& observation_frequency) {
  vec rho_cluster = rho_old.col(cluster_index);
  rho_old.col(cluster_index) = make_new_rho(rho_cluster, rankings, alpha_old(cluster_index),
              leap_size, metric, observation_frequency);

  // Save rho if appropriate
  if(t % rho_thinning == 0){
    if(cluster_index == 0) ++rho_index;
    rho.slice(rho_index).col(cluster_index) = rho_old.col(cluster_index);
  }
}

void parameters::update_shape(int t, const mat& rankings,
                              const Rcpp::List& constraints) {


  int n_assessors = rankings.n_cols;
  int sum_1{};
  int sum_2{};
  for(int i = 0; i < n_assessors; ++i){
    Rcpp::List assessor_constraints = Rcpp::as<Rcpp::List>(constraints[i]);
    for(int j = 0; j < n_items; ++j) {
      uvec items_above = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[j]);
      uvec items_below = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[j]);

      for(unsigned int k = 0; k < items_above.n_elem; ++k){
        int g = (as_scalar(rankings.col(i).row(j)) < as_scalar(rankings.col(i).row(items_above(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
      for(unsigned int k = 0; k < items_below.n_elem; ++k){
        int g = (as_scalar(rankings.col(i).row(j)) > as_scalar(rankings.col(i).row(items_below(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
    }
  }
  shape_1(t) = kappa_1 + sum_1;
  shape_2(t) = kappa_2 + sum_2;
  theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
}

void parameters::update_alpha(
    int cluster_index,
    int alpha_index,
    const mat& rankings,
    const vec& observation_frequency,
    const Rcpp::List& logz_list) {

  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                                   std::log(alpha_old(cluster_index)));

  double rank_dist = rank_dist_sum(rankings, rho_old.col(cluster_index), metric, observation_frequency);

  // Difference between current and proposed alpha
  double alpha_diff = alpha_old(cluster_index) - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    lambda * alpha_diff +
    sum(observation_frequency) * (
        get_partition_function(n_items, alpha_old(cluster_index), logz_list, metric) -
          get_partition_function(n_items, alpha_proposal, logz_list, metric)
    ) + std::log(alpha_proposal) - std::log(alpha_old(cluster_index));

  // Draw a uniform random number
  double u = std::log(randu<double>());

  if(ratio > u){
    alpha(cluster_index, alpha_index) = alpha_proposal;
  } else {
    alpha(cluster_index, alpha_index) = alpha_old(cluster_index);
  }
}
