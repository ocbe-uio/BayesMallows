#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void smc_mallows_new_users_augment_partial(
    arma::cube& aug_rankings,
    arma::vec& aug_prob,
    const arma::mat& rho_samples,
    const arma::vec& alpha_samples,
    const int& num_obs,
    const int& num_new_obs,
    const arma::mat& rankings,
    const std::string& aug_method,
    const double& alpha,
    const bool& augment_alpha,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_cols;
  int n_items = rho_samples.n_rows;
  ivec ranks = regspace<ivec>(1, n_items);
  for (int ii{}; ii < n_particles; ++ii) {
    for (int jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
      vec partial_ranking = rankings.row(jj).t();

      // find items missing from original observed ranking
      const uvec& unranked_items = find_nonfinite(partial_ranking);

      // find ranks missing from ranking
      const vec& missing_ranks = setdiff_template(ranks, partial_ranking);

      // fill in missing ranks based on choice of augmentation method
      Rcpp::List proposal{};
      const bool pseudo = is_pseudo(aug_method, metric);
      if (!pseudo) {
        // create new augmented ranking by sampling remaining ranks from set uniformly
        partial_ranking.elem(find_nonfinite(partial_ranking)) = shuffle(missing_ranks);

        aug_rankings(span(jj), span::all, span(ii)) = partial_ranking;
        aug_prob(ii) = divide_by_fact(aug_prob(ii), missing_ranks.n_elem);
      } else {
        // randomly permute the unranked items to give the order in which they will be allocated
        uvec item_ordering = shuffle(unranked_items);
        const vec& rho_s = rho_samples.col(ii);
        proposal = calculate_forward_probability(
          item_ordering, partial_ranking, missing_ranks, rho_s,
          augment_alpha ? alpha_samples(ii) : alpha, n_items, metric
        );

        const vec& a_rank = proposal["aug_ranking"];
        const double& f_prob = proposal["forward_prob"];
        aug_rankings(span(jj), span::all, span(ii)) = a_rank;
        aug_prob(ii) = aug_prob(ii) * f_prob;
      }
    }
  }
}

vec initialize_alpha(
    const int& n_particles,
    const Rcpp::Nullable<vec>& alpha_init){
  if(alpha_init.isNotNull()) {
    return Rcpp::as<vec>(alpha_init);
  } else {
    return Rcpp::rexp(n_particles, 1);
  }

}

void smc_mallows_new_users_resample(
    mat& rho_samples, vec& alpha_samples, cube& aug_rankings,
    const vec& norm_wgt, const int& num_obs,
    const bool& augment_alpha, const bool& partial
){
  int n_particles = rho_samples.n_cols;
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);

  rho_samples = rho_samples.cols(index);

  /* Replacing tt + 1 column on alpha_samples ------------- */
  if(augment_alpha){
    alpha_samples = alpha_samples.rows(index);
  }

  if(partial){
    cube aug_rankings_index = aug_rankings.slices(index);
    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(span(0, num_obs - 1), span::all, span::all);
  }
}

void smc_mallows_new_users_reweight(
    vec& log_inc_wgt,
    double& effective_sample_size,
    vec& norm_wgt,
    const cube& aug_rankings,
    const mat& observed_rankings,
    const mat& rho_samples,
    const double& alpha,
    const vec& alpha_samples,
    const Rcpp::Nullable<vec> logz_estimate,
    const Rcpp::Nullable<vec> cardinalities,
    const int& num_obs,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& augment_alpha,
    const bool& partial,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_cols;
  int n_items = rho_samples.n_rows;
  for (int ii{}; ii < n_particles; ++ii) {


    double alpha_used = augment_alpha ? alpha_samples(ii) : alpha;


    /* Calculating log_z_alpha and log_likelihood ----------- */
    const double log_z_alpha = get_partition_function(
      n_items, alpha_used, cardinalities, logz_estimate, metric
    );

    double log_likelihood = get_exponent_sum(                          \
      alpha_used, rho_samples.col(ii), n_items,
      partial ? aug_rankings(span(num_obs - num_new_obs, num_obs - 1), span::all, span(ii)) : observed_rankings,
      metric
    );
    log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
  }

  /* normalise weights ------------------------------------ */
  norm_wgt = normalize_weights(log_inc_wgt);

  effective_sample_size = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);
}
