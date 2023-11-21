#include "parameters.h"

using namespace arma;

Data::Data(
  const Rcpp::List& data,
  const Rcpp::List& compute_options
) :
  rankings { Rcpp::as<mat>(data["rankings"]).t() },
  constraints { Rcpp::as<Rcpp::List>(data["constraints"]) },
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  augpair { constraints.length() > 0 },
  any_missing { !is_finite(rankings) },
  save_aug { Rcpp::as<bool>(compute_options["save_aug"]) },
  aug_thinning { Rcpp::as<unsigned int>(compute_options["aug_thinning"]) }
  {

    if(any_missing){
      set_up_missing(rankings, missing_indicator);
      initialize_missing_ranks(rankings, missing_indicator);
    }

    if(save_aug){
      unsigned int nmc = Rcpp::as<unsigned int>(compute_options["nmc"]);
      augmented_data.set_size(n_items, n_assessors, std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
      augmented_data.slice(0) = rankings;
    }

  }

Priors::Priors(
  const Rcpp::List& priors
) : lambda { verify_positive(Rcpp::as<double>(priors["lambda"])) },
  kappa_1 { Rcpp::as<unsigned int>(priors["kappa_1"]) },
  kappa_2 { Rcpp::as<unsigned int>(priors["kappa_2"]) }
  {

}

Parameters::Parameters(
  const Rcpp::List& model,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values,
  const unsigned int n_items) :
  alpha_jump { Rcpp::as<int>(compute_options["alpha_jump"]) },
  alpha_prop_sd { verify_positive(Rcpp::as<double>(compute_options["alpha_prop_sd"])) },
  error_model { verify_error_model(Rcpp::as<std::string>(model["error_model"])) },
  leap_size { Rcpp::as<int>(compute_options["leap_size"]) },
  metric { verify_metric(Rcpp::as<std::string>(model["metric"])) },
  n_clusters { Rcpp::as<int>(model["n_clusters"]) },
  nmc { Rcpp::as<int>(compute_options["nmc"]) },
  rho_thinning { Rcpp::as<int>(compute_options["rho_thinning"]) }
  {

    alpha.set_size(n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / alpha_jump)));
    double alpha_init = initial_values["alpha_init"];
    alpha.col(0).fill(alpha_init);
    alpha_old = alpha.col(0);

    rho.set_size(n_items, n_clusters, std::ceil(static_cast<double>(nmc * 1.0 / rho_thinning)));
    Rcpp::Nullable<mat> rho_init = initial_values["rho_init"];
    rho.slice(0) = initialize_rho(n_items, n_clusters, rho_init);
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

void Parameters::update_rho(int cluster_index, int t, int& rho_index,
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

void Parameters::update_shape(int t, const mat& rankings,
                              const Rcpp::List& constraints,
                              const Priors& priors) {

  const unsigned int n_items = rankings.n_rows;
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

  shape_1(t) = priors.kappa_1 + sum_1;
  shape_2(t) = priors.kappa_2 + sum_2;
  theta(t) = rtruncbeta(shape_1(t), shape_2(t), 0.5);
}

void Parameters::update_alpha(
    int cluster_index,
    int alpha_index,
    const mat& rankings,
    const vec& observation_frequency,
    const Rcpp::List& logz_list,
    const Priors& priors) {

  const unsigned int n_items = rankings.n_rows;
  double alpha_proposal = std::exp(randn<double>() * alpha_prop_sd +
                                   std::log(alpha_old(cluster_index)));

  double rank_dist = rank_dist_sum(rankings, rho_old.col(cluster_index), metric, observation_frequency);

  // Difference between current and proposed alpha
  double alpha_diff = alpha_old(cluster_index) - alpha_proposal;

  // Compute the Metropolis-Hastings ratio
  double ratio =
    alpha_diff / n_items * rank_dist +
    priors.lambda * alpha_diff +
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
