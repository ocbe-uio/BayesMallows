#include "classes.h"
#include "missing_data.h"
#include "parameterupdates.h"
#include "misc.h"
#include "partitionfuns.h"
#include "sample.h"
#include "distances.h"
#include "pairwise_comparisons.h"

using namespace arma;

static std::string verify_metric(const std::string input) {
  bool check = (input.compare("footrule") == 0) ||
    (input.compare("spearman") == 0) ||
    (input.compare("cayley") == 0) ||
    (input.compare("kendall") == 0) ||
    (input.compare("ulam") == 0) ||
    (input.compare("hamming") == 0);
  if(!check) {
    Rcpp::stop("Unknown metric.\n");
  }
  return input;
}

static std::string verify_error_model(const std::string input) {
  bool check = (input.compare("none") == 0) ||
    (input.compare("bernoulli") == 0);
  if(!check) {
    Rcpp::stop("Unknown error model.\n");
  }
  return input;
}

Data::Data(
  const Rcpp::List& data
) :
  rankings { Rcpp::as<mat>(data["rankings"]).t() },
  constraints { Rcpp::as<Rcpp::List>(data["constraints"]) },
  n_assessors { rankings.n_cols },
  n_items { rankings.n_rows },
  observation_frequency { Rcpp::as<vec>(data["observation_frequency"]) } {}

SMCData::SMCData(
  const Rcpp::List& data,
  const Rcpp::List& new_data
) : Data(data),
  new_rankings { Rcpp::as<mat>(new_data["rankings"]).t() },
  num_new_obs { new_rankings.n_cols },
  consistent {Rcpp::as<umat>(new_data["consistent"])} {}

SMCParameters::SMCParameters(
  const Rcpp::List& model_options,
  const Rcpp::List& smc_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values
) :
  n_particles { Rcpp::as<unsigned int>(smc_options["n_particles"]) },
  mcmc_steps { Rcpp::as<unsigned int>(smc_options["mcmc_steps"]) },
  alpha_samples { Rcpp::as<arma::vec>(initial_values["alpha_init"]) },
  rho_samples { Rcpp::as<arma::mat>(initial_values["rho_init"]) },
  alpha_prop_sd { verify_positive(Rcpp::as<double>(compute_options["alpha_prop_sd"])) },
  leap_size { Rcpp::as<unsigned int>(compute_options["leap_size"]) },
  metric { verify_metric(Rcpp::as<std::string>(model_options["metric"])) },
  norm_wgt { ones(n_particles) } {}

Augmentation::Augmentation(
  Data& dat,
  const Rcpp::List& compute_options
) :
  augpair { dat.constraints.length() > 0 },
  any_missing { !is_finite(dat.rankings) },
  save_aug { Rcpp::as<bool>(compute_options["save_aug"]) },
  aug_thinning { Rcpp::as<unsigned int>(compute_options["aug_thinning"]) },
  swap_leap { Rcpp::as<unsigned int>(compute_options["swap_leap"]) } {

  if(any_missing){
    set_up_missing(dat.rankings, missing_indicator);
    initialize_missing_ranks(dat.rankings, missing_indicator);
  }

  if(save_aug){
    unsigned int nmc = Rcpp::as<unsigned int>(compute_options["nmc"]);
    augmented_data.set_size(dat.n_items, dat.n_assessors,
                            std::ceil(static_cast<double>(nmc * 1.0 / aug_thinning)));
    augmented_data.slice(0) = dat.rankings;
  }}

SMCAugmentation::SMCAugmentation(
  SMCData& dat,
  const Rcpp::List& smc_options,
  const Rcpp::List& initial_values,
  const unsigned int n_particles
) :
  aug_method { Rcpp::as<std::string>(smc_options["aug_method"]) },
  aug_prob { arma::ones(n_particles) },
  any_missing { !is_finite(dat.rankings) },
  aug_init { Rcpp::as<Rcpp::Nullable<arma::cube>>(initial_values["aug_init"])}
  {
    if(any_missing){
      set_up_missing(dat.rankings, missing_indicator);
      augmented_data.set_size(dat.n_items, dat.n_assessors, n_particles);

      for(int i{}; i < n_particles; i++) {
        augmented_data.slice(i) = dat.rankings;
        initialize_missing_ranks(augmented_data.slice(i), missing_indicator);
      }

      if(aug_init.isNotNull()) {
        augmented_data(
          span::all,
          span(0, dat.rankings.n_cols - dat.new_rankings.n_cols - 1),
          span::all) = Rcpp::as<cube>(aug_init);
      }
    }
  }

Priors::Priors(
  const Rcpp::List& priors
) : lambda { verify_positive(Rcpp::as<double>(priors["lambda"])) },
  kappa_1 { Rcpp::as<unsigned int>(priors["kappa_1"]) },
  kappa_2 { Rcpp::as<unsigned int>(priors["kappa_2"]) },
  psi { Rcpp::as<unsigned int>(priors["psi"]) } {}

Parameters::Parameters(
  const Rcpp::List& model_options,
  const Rcpp::List& compute_options,
  const Rcpp::List& initial_values,
  const unsigned int n_items) :
  n_clusters { Rcpp::as<int>(model_options["n_clusters"]) },
  nmc { Rcpp::as<int>(compute_options["nmc"]) },
  metric { verify_metric(Rcpp::as<std::string>(model_options["metric"])) },
  error_model { verify_error_model(Rcpp::as<std::string>(model_options["error_model"])) },
  alpha_jump { Rcpp::as<int>(compute_options["alpha_jump"]) },
  alpha_prop_sd { verify_positive(Rcpp::as<double>(compute_options["alpha_prop_sd"])) },
  element_indices { regspace<uvec>(0, n_items - 1) },
  leap_size { Rcpp::as<int>(compute_options["leap_size"]) },
  rho_thinning { Rcpp::as<int>(compute_options["rho_thinning"]) }
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
      for (int i = 0; i < n_clusters; ++i) {
        rho.slice(0).col(i) = randperm<vec>(n_items) + 1;
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

Clustering::Clustering(const Parameters& pars,
                       const Rcpp::List& compute_options,
                       const unsigned int n_assessors) :
  index { regspace<uvec>(0, pars.n_clusters - 1) },
  clustering {pars.n_clusters > 1},
  clus_thinning { Rcpp::as<unsigned int>(compute_options["clus_thinning"]) },
  include_wcd { Rcpp::as<bool>(compute_options["include_wcd"]) },
  save_ind_clus { Rcpp::as<bool>(compute_options["save_ind_clus"]) } {
    int n_cluster_assignments = pars.n_clusters > 1 ? std::ceil(static_cast<double>(pars.nmc * 1.0 / clus_thinning)) : 1;
    cluster_probs.set_size(pars.n_clusters, n_cluster_assignments);
    cluster_probs.col(0).fill(1.0 / pars.n_clusters);
    current_cluster_probs = cluster_probs.col(0);
    cluster_assignment.set_size(n_assessors, n_cluster_assignments);
    cluster_assignment.col(0) = randi<uvec>(n_assessors, distr_param(0, pars.n_clusters - 1));
    current_cluster_assignment = cluster_assignment.col(0);

    dist_mat.set_size(n_assessors, pars.n_clusters);

    within_cluster_distance.set_size(pars.n_clusters, include_wcd ? pars.nmc : 1);
    update_wcd(0);
}

void Parameters::update_rho(int t, int& rho_index, const Data& dat,
                            const uvec& current_cluster_assignment) {
  for(int i{}; i < n_clusters; ++i){
    const uvec cluster_indicator = find(current_cluster_assignment == i);
    const mat cluster_rankings = dat.rankings.submat(
      element_indices, cluster_indicator);
    const vec cluster_frequency =
      dat.observation_frequency.elem(cluster_indicator);

    vec rho_cluster = rho_old.col(i);
    rho_old.col(i) = make_new_rho(rho_cluster, cluster_rankings, alpha_old(i),
                leap_size, metric, cluster_frequency);

    if(t % rho_thinning == 0){
      if(i == 0) ++rho_index;
      rho.slice(rho_index).col(i) = rho_old.col(i);
    }
  }
}

void SMCParameters::update_rho(
    const unsigned int particle_index,
    const SMCData& dat) {
  rho_samples.col(particle_index) = make_new_rho(
    rho_samples.col(particle_index), dat.rankings,
    alpha_samples(particle_index), leap_size, metric,
    dat.observation_frequency);
}

void Parameters::update_shape(int t, const Data& dat,
                              const Priors& priors) {
  if(error_model != "bernoulli") return;
  int sum_1{};
  int sum_2{};
  for(int i = 0; i < dat.n_assessors; ++i){
    Rcpp::List assessor_constraints = Rcpp::as<Rcpp::List>(dat.constraints[i]);
    for(int j = 0; j < dat.n_items; ++j) {
      uvec items_above = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[j]);
      uvec items_below = Rcpp::as<uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[j]);

      for(unsigned int k = 0; k < items_above.n_elem; ++k){
        int g = (as_scalar(dat.rankings.col(i).row(j)) < as_scalar(dat.rankings.col(i).row(items_above(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
      for(unsigned int k = 0; k < items_below.n_elem; ++k){
        int g = (as_scalar(dat.rankings.col(i).row(j)) > as_scalar(dat.rankings.col(i).row(items_below(k) - 1)));
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
    int alpha_index,
    const Data& dat,
    const Rcpp::List& logz_list,
    const Priors& priors,
    const uvec& current_cluster_assignment) {

  for(int i = 0; i < n_clusters; ++i) {
    const uvec cluster_indicator = find(current_cluster_assignment == i);
    const mat cluster_rankings = dat.rankings.submat(
      element_indices, cluster_indicator);
    const vec cluster_frequency =
      dat.observation_frequency.elem(cluster_indicator);

    AlphaRatio test = make_new_alpha(
      alpha_old(i), rho_old.col(i),
      alpha_prop_sd, metric, logz_list, cluster_rankings,
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
    const Rcpp::List& logz_list,
    const Priors& priors) {

  AlphaRatio test = make_new_alpha(
    alpha_samples(particle_index),
    rho_samples.col(particle_index),
    alpha_prop_sd, metric, logz_list,
    dat.rankings, dat.observation_frequency, dat.n_items, priors
  );

  if(test.accept){
    alpha_samples(particle_index) = test.proposal;
  }
}

void Clustering::update_cluster_probs(const Parameters& pars, const Priors& pris){
  vec cluster_probs(pars.n_clusters);

  for(int i = 0; i < pars.n_clusters; ++i){
    cluster_probs(i) = R::rgamma(sum(current_cluster_assignment == i) + pris.psi, 1.0);
  }
  current_cluster_probs = normalise(cluster_probs, 1);
}

void Clustering::update_cluster_labels(
    const int t,
    const Data& dat,
    const Parameters& pars,
    const Rcpp::List& logz_list
){
  uvec new_cluster_assignment(dat.n_assessors);
  mat assignment_prob(dat.n_assessors, pars.n_clusters);
  for(int i = 0; i < pars.n_clusters; ++i){
    // Compute the logarithm of the unnormalized probability
    assignment_prob.col(i) = std::log(cluster_probs(i)) -
      pars.alpha_old(i) / dat.n_items * dist_mat.col(i) -
      get_partition_function(dat.n_items, pars.alpha_old(i), logz_list, pars.metric);
  }

  for(int i = 0; i < dat.n_assessors; ++i){
    // Exponentiate to get unnormalized prob relative to max
    rowvec probs = exp(assignment_prob.row(i) -
      max(assignment_prob.row(i)));

    // Normalize with 1-norm
    assignment_prob.row(i) = normalise(probs, 1);
    new_cluster_assignment(span(i)) = sample(index, 1, false, assignment_prob.row(i).t());
  }

  if(save_ind_clus){
    assignment_prob.save(std::string("cluster_probs") + std::to_string(t + 1) + std::string(".csv"), csv_ascii);
  }
  current_cluster_assignment = new_cluster_assignment;
}

void Clustering::update_wcd(const int t){
  if(!include_wcd) return;

  int n_clusters = dist_mat.n_cols;
  vec wcd(n_clusters);

  for(int i = 0; i < n_clusters; ++i){
    mat dist_vec = dist_mat.submat(find(current_cluster_assignment == i), index.subvec(i, i));
    wcd(i) = accu(dist_vec);
  }
  within_cluster_distance.col(t) = wcd;
}

void Clustering::update_dist_mat(const Data& dat, const Parameters& pars){
  if(clustering | include_wcd) {
    for(int i = 0; i < pars.n_clusters; ++i)
      dist_mat.col(i) = rank_dist_vec(dat.rankings, pars.rho_old.col(i),
                   pars.metric, dat.observation_frequency);
  }
}

void Augmentation::augment_pairwise(
    const unsigned int t,
    Data& dat,
    const Parameters& pars,
    const Clustering& clus
){
  if(!augpair) return;
  for(int i = 0; i < dat.n_assessors; ++i) {
    vec proposal;
    int g_diff = 0;

    if(pars.error_model == "none"){
      proposal = propose_pairwise_augmentation(dat.rankings.col(i), Rcpp::as<Rcpp::List>(dat.constraints[i]));
    } else if(pars.error_model == "bernoulli"){
      proposal = propose_swap(dat.rankings.col(i), Rcpp::as<Rcpp::List>(dat.constraints[i]), g_diff, swap_leap);
    } else {
      Rcpp::stop("error_model must be 'none' or 'bernoulli'");
    }

    double u = std::log(randu<double>());

    // Find which cluster the assessor belongs to
    int cluster = clus.current_cluster_assignment(i);

    double ratio = -pars.alpha_old(cluster) / dat.n_items *
      (get_rank_distance(proposal, pars.rho_old.col(cluster), pars.metric) -
      get_rank_distance(dat.rankings.col(i), pars.rho_old.col(cluster), pars.metric));

    if(pars.error_model != "none") {
      ratio += g_diff * std::log(pars.theta(t) / (1 - pars.theta(t)));
    }

    if(ratio > u){
      dat.rankings.col(i) = proposal;
    }
  }
}

void Augmentation::update_missing_ranks(
    Data& dat,
    const Clustering& clus,
    const Parameters& pars) {
  if(!any_missing) return;

  for(int i = 0; i < dat.n_assessors; ++i){
    int cluster = clus.current_cluster_assignment(i);

    dat.rankings.col(i) = make_new_augmentation(
      dat.rankings.col(i), missing_indicator.col(i),
      pars.alpha_old(cluster), pars.rho_old.col(cluster),
      pars.metric
    );
  }
}

void SMCAugmentation::reweight(
  SMCParameters& pars,
  const SMCData& dat,
  const Rcpp::List& logz_list
) {
  augment_partial(pars, dat);

  vec log_inc_wgt(pars.n_particles);
  for (size_t particle{}; particle < pars.n_particles; ++particle) {

    const double log_z_alpha = get_partition_function(
      dat.n_items, pars.alpha_samples(particle), logz_list, pars.metric
    );
    double item_correction_contribution{};
    if(!dat.consistent.is_empty()) {
      for(size_t user{}; user < dat.n_assessors - dat.num_new_obs; user++) {
        if(dat.consistent(user, particle) == 0) {
          vec previous_augmented_ranking = augmented_data(span::all, span(user), span(particle - 1));
          vec current_augmented_ranking = augmented_data(span::all, span(user), span(particle));

          double previous_distance = get_rank_distance(
            previous_augmented_ranking, pars.rho_samples.col(particle),
            pars.metric);
          double current_distance = get_rank_distance(
            current_augmented_ranking, pars.rho_samples.col(particle),
            pars.metric);

          item_correction_contribution -= pars.alpha_samples(particle) / dat.n_items *
            (current_distance - previous_distance);

        }
      }
    }

    double new_user_contribution{};
    if(dat.num_new_obs > 0) {
      const mat new_rankings = !any_missing ? dat.new_rankings :
      augmented_data(
        span::all,
        span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1),
        span(particle));

      new_user_contribution = -pars.alpha_samples(particle) / dat.n_items *
        rank_dist_sum(new_rankings, pars.rho_samples.col(particle), pars.metric,
                      dat.observation_frequency(span(dat.n_assessors - dat.num_new_obs, dat.n_assessors - 1)));
    }

    log_inc_wgt(particle) = new_user_contribution + item_correction_contribution -
      dat.num_new_obs * log_z_alpha - log(aug_prob(particle));
  }

  pars.norm_wgt = normalize_weights(log_inc_wgt);
  pars.effective_sample_size = 1 / std::pow(norm(pars.norm_wgt, 2), 2);
}

void SMCAugmentation::augment_partial(
    const SMCParameters& pars,
    const SMCData& dat
){
  if(!any_missing) return;
  for (size_t particle{}; particle < pars.n_particles; ++particle) {

    for (size_t user{}; user < dat.n_assessors; ++user) {
      if(user < dat.n_assessors - dat.num_new_obs) {
        if(dat.consistent.is_empty()) continue;
        if(dat.consistent(user, particle) == 1) continue;
      }

      uvec unranked_items = shuffle(find(missing_indicator.col(user) == 1));

      if (aug_method != "pseudo") {
        augmented_data(span::all, span(user), span(particle)) =
          propose_augmentation(augmented_data(span::all, span(user), span(particle)),
                               missing_indicator.col(user));

      } else {
        PseudoProposal pprop = make_pseudo_proposal(
          unranked_items, augmented_data(span::all, span(user), span(particle)),
          pars.alpha_samples(particle), pars.rho_samples.col(particle), pars.metric
        );

        augmented_data(span::all, span(user), span(particle)) = pprop.rankings;
        aug_prob(particle) *= pprop.probability;
      }
    }
  }
}

void SMCAugmentation::update_data(
    const unsigned int particle_index, SMCData& dat) {
  if(!any_missing) return;
  dat.rankings = augmented_data.slice(particle_index);
}

uvec SMCParameters::draw_resampling_index() {
  return sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);
}

void SMCParameters::resample(const uvec& index) {
  rho_samples = rho_samples.cols(index);
  alpha_samples = alpha_samples.rows(index);
}

void SMCAugmentation::resample(const uvec& index) {
  if(!any_missing) return;
  augmented_data = augmented_data.slices(index);
}

void SMCAugmentation::update_missing_ranks(
    const unsigned int particle_index,
    const SMCData& dat,
    const SMCParameters& pars) {
  if(!any_missing) return;

  for (int jj = dat.n_assessors - dat.num_new_obs; jj < dat.n_assessors; ++jj) {
    augmented_data(span::all, span(jj), span(particle_index)) = make_new_augmentation(
      augmented_data(span::all, span(jj), span(particle_index)),
      missing_indicator.col(jj),
      pars.alpha_samples(particle_index),
      pars.rho_samples.col(particle_index),
      pars.metric, aug_method == "pseudo"
    );
  }
}
