#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"
#include "missing_data.h"
#include "distances.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

Rcpp::List make_pseudo_proposal(
    uvec unranked_items, vec rankings, const double& alpha, const vec& rho,
    const std::string metric, const bool forward = true
) {

  int n_items = rankings.n_elem;

  double prob = 1;
  while(unranked_items.n_elem > 0) {
    vec available_rankings = rankings(unranked_items);
    int item_to_rank = unranked_items(0);

    vec log_numerator(available_rankings.n_elem);
    for(int ll{}; ll < available_rankings.n_elem; ll++) {
      log_numerator(ll) = -alpha / n_items *
        get_rank_distance(rho(span(item_to_rank)), available_rankings(span(ll)), metric);
    }
    vec sample_probs = normalize_weights(log_numerator);
    if(forward) {
      rankings(span(item_to_rank)) =
        sample(available_rankings, 1, false, sample_probs);
    }


    int ranking_chosen = as_scalar(find(rankings(item_to_rank) == available_rankings));

    prob *= sample_probs(ranking_chosen);
    if(available_rankings.n_elem <= 1) break;
    unranked_items = unranked_items.subvec(1, available_rankings.n_elem - 1);

    rankings(unranked_items) = setdiff_template(available_rankings, available_rankings(span(ranking_chosen)));
  }

  return Rcpp::List::create(
    Rcpp::Named("proposal") = rankings,
    Rcpp::Named("probability") = prob
  );
}

void smc_mallows_new_users_augment_partial(
    arma::cube& augmented_data,
    arma::vec& aug_prob,
    const arma::mat& rho_samples,
    const arma::vec& alpha_samples,
    const int& num_new_obs,
    const std::string& aug_method,
    const umat& missing_indicator,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_cols;
  int num_obs = augmented_data.n_cols;

  for (int ii{}; ii < n_particles; ++ii) {
    double alpha = alpha_samples(ii);
    vec rho = rho_samples.col(ii);
    for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
      uvec unranked_items = shuffle(find(missing_indicator.col(jj) == 1));

      if (aug_method != "pseudo") {
        augmented_data(span::all, span(jj), span(ii)) =
          propose_augmentation(augmented_data(span::all, span(jj), span(ii)),
                               missing_indicator.col(jj));

        aug_prob(ii) = divide_by_fact(aug_prob(ii), unranked_items.n_elem);

      } else {
        Rcpp::List pprop = make_pseudo_proposal(
          unranked_items, augmented_data(span::all, span(jj), span(ii)),
          alpha, rho, metric
        );

        vec ar = pprop["proposal"];
        augmented_data(span::all, span(jj), span(ii)) = ar;
        double prob = pprop["probability"];
        aug_prob(ii) *= prob;
      }
    }
  }
}

void smc_mallows_new_users_resample(
    mat& rho_samples, vec& alpha_samples, cube& augmented_data,
    const vec& norm_wgt,
    const bool& partial
){
  int n_particles = rho_samples.n_cols;
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);
  rho_samples = rho_samples.cols(index);
  alpha_samples = alpha_samples.rows(index);

  if(partial) augmented_data = augmented_data.slices(index);
}

void smc_mallows_new_users_reweight(
    vec& log_inc_wgt,
    double& effective_sample_size,
    vec& norm_wgt,
    const cube& augmented_data,
    const mat& observed_rankings,
    const mat& rho_samples,
    const vec& alpha_samples,
    const Rcpp::Nullable<vec> logz_estimate,
    const Rcpp::Nullable<vec> cardinalities,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& partial,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_cols;
  int n_items = rho_samples.n_rows;
  int num_obs = augmented_data.n_cols;
  for (int ii{}; ii < n_particles; ++ii) {

    const double log_z_alpha = get_partition_function(
      n_items, alpha_samples(ii), cardinalities, logz_estimate, metric
    );

    mat rankings = !partial ? observed_rankings : augmented_data(span::all,
      span(num_obs - num_new_obs, num_obs - 1),
      span(ii));

    double log_likelihood = -alpha_samples(ii) / n_items *
      rank_dist_sum(rankings, rho_samples.col(ii), metric, ones(rankings.n_cols));

    log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
  }

  norm_wgt = normalize_weights(log_inc_wgt);
  effective_sample_size = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);
}


