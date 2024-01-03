#include "distributions.h"
#include "classes.h"
#include "distances.h"
#include "leapandshift.h"
#include "missing_data.h"
using namespace arma;

struct AlphaRatio{
  AlphaRatio(double proposal, bool accept) :
  proposal { proposal }, accept { accept } {}
  ~AlphaRatio() = default;
  double proposal;
  bool accept;
};

AlphaRatio make_new_alpha(
    const double& alpha_old,
    const vec& rho_old,
    const double& alpha_prop_sd,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartitionFunction>& pfun,
    const arma::mat& rankings,
    const arma::vec& observation_frequency,
    const double& n_items,
    const Priors& priors) {

  double alpha_proposal = R::rlnorm(std::log(alpha_old), alpha_prop_sd);
  double rank_dist = sum(distfun->d(rankings, rho_old) % observation_frequency);
  double alpha_diff = alpha_old - alpha_proposal;

  double ratio =
    alpha_diff / n_items * rank_dist +
    priors.lambda * alpha_diff +
    sum(observation_frequency) *
    (pfun->logz(alpha_old) - pfun->logz(alpha_proposal)) +
    std::log(alpha_proposal) - std::log(alpha_old);

  return AlphaRatio{alpha_proposal, ratio > std::log(R::unif_rand())};
}

vec make_new_rho(
    vec current_rho,
    const mat& rankings,
    double alpha_old,
    int leap_size,
    const std::unique_ptr<Distance>& distfun,
    vec observation_frequency) {

  vec rho_proposal;
  uvec indices;
  double prob_backward, prob_forward;
  int n_items = current_rho.n_elem;

  leap_and_shift(
    rho_proposal, indices, prob_backward, prob_forward,
    current_rho, leap_size, distfun);

  const arma::mat& r = rankings.rows(indices);
  double dist_new = arma::sum(
    distfun->d(r, rho_proposal(indices)) % observation_frequency
  );
  double dist_old = arma::sum(
    distfun->d(r, current_rho(indices)) % observation_frequency
  );

  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  if(ratio > std::log(R::unif_rand())){
    return(rho_proposal);
  } else {
    return(current_rho);
  }
}


Parameters::Parameters(
  const Rcpp::List& model_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values,
  const unsigned int n_items) :
  n_clusters { model_options["n_clusters"] },
  nmc { compute_options["nmc"] },
  error_model(model_options["error_model"]),
  alpha_jump { compute_options["alpha_jump"] },
  element_indices { regspace<uvec>(0, n_items - 1) },
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  leap_size { compute_options["leap_size"] },
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
        Rcpp::IntegerVector a = Rcpp::sample(n_items, n_items);
        rho.slice(0).col(i) = Rcpp::as<vec>(Rcpp::wrap(a));
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

SMCParameters::SMCParameters(
  const Rcpp::List& model_options,
  const Rcpp::List& smc_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values
) :
  n_particles { smc_options["n_particles"] },
  mcmc_steps { smc_options["mcmc_steps"] },
  alpha_samples(initial_values["alpha_init"]) ,
  rho_samples(initial_values["rho_init"]),
  alpha_prop_sd { compute_options["alpha_prop_sd"] },
  leap_size { compute_options["leap_size"] },
  log_inc_wgt { zeros(n_particles) } {}

void Parameters::update_rho(
    int t,
    int& rho_index,
    const Data& dat,
    const uvec& current_cluster_assignment,
    const std::unique_ptr<Distance>& distfun
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

void SMCParameters::update_rho(
    const unsigned int particle_index,
    const SMCData& dat,
    const std::unique_ptr<Distance>& distfun) {
  rho_samples.col(particle_index) = make_new_rho(
    rho_samples.col(particle_index), dat.rankings,
    alpha_samples(particle_index), leap_size, distfun,
    dat.observation_frequency);
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

void SMCParameters::update_alpha(
    const unsigned int particle_index,
    const SMCData& dat,
    const std::unique_ptr<PartitionFunction>& pfun,
    const std::unique_ptr<Distance>& distfun,
    const Priors& priors) {

  AlphaRatio test = make_new_alpha(
    alpha_samples(particle_index),
    rho_samples.col(particle_index),
    alpha_prop_sd, distfun, pfun,
    dat.rankings, dat.observation_frequency,
    dat.n_items, priors
  );
  if(test.accept){
    alpha_samples(particle_index) = test.proposal;
  }
}

uvec SMCParameters::draw_resampling_index() {
  Rcpp::IntegerVector inds = Rcpp::seq(0, log_inc_wgt.size() - 1);
  vec norm_wgt = exp(log_inc_wgt - max(log_inc_wgt) -
    log(sum(exp(log_inc_wgt - max(log_inc_wgt)))));
  Rcpp::NumericVector probs = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(norm_wgt));
  Rcpp::IntegerVector result = Rcpp::sample(inds, log_inc_wgt.size(), true, probs);
  return Rcpp::as<arma::uvec>(result);
}

void SMCParameters::resample(const uvec& index) {
  rho_samples = rho_samples.cols(index);
  alpha_samples = alpha_samples.rows(index);
}

void SMCAugmentation::resample(const uvec& index, const SMCData& dat) {
  if(!dat.any_missing) return;
  augmented_data = augmented_data.slices(index);
}

void SMCAugmentation::update_missing_ranks(
    const unsigned int particle_index,
    const SMCData& dat,
    const SMCParameters& pars,
    const std::unique_ptr<Distance>& distfun) {
  if(!dat.any_missing) return;

  for (unsigned int jj = dat.n_assessors - dat.num_new_obs;
       jj < dat.n_assessors; ++jj) {
    augmented_data(span::all, span(jj), span(particle_index)) =
      make_new_augmentation(
        augmented_data(span::all, span(jj), span(particle_index)),
        missing_indicator.col(jj),
        pars.alpha_samples(particle_index),
        pars.rho_samples.col(particle_index),
        distfun, log_aug_prob(jj, particle_index), aug_method == "pseudo"
    );
  }
}
