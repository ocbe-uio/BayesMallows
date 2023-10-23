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
    const arma::cube& rho_samples,
    const arma::mat& alpha_samples,
    const int& num_obs,
    const int& num_new_obs,
    const arma::mat& rankings,
    const std::string& aug_method,
    const int& tt,
    const double& alpha,
    const bool& augment_alpha,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_rows;
  int n_items = rho_samples.n_cols;
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
        const rowvec& rho_s = rho_samples(span(ii), span::all, span(tt + 1));
        proposal = calculate_forward_probability(
          item_ordering, partial_ranking, missing_ranks, rho_s.t(),
          augment_alpha ? alpha_samples(ii, tt + 1) : alpha, n_items, metric
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
    cube& rho_samples, mat& alpha_samples, cube& aug_rankings,
    const vec& norm_wgt, const int& tt, const int& num_obs,
    const bool& augment_alpha, const bool& partial
){
  int n_particles = rho_samples.n_rows;
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, n_particles - 1), n_particles, true, norm_wgt);
  uvec tt_vec;
  tt_vec = tt;

  /* Replacing tt + 1 slice on rho_samples ---------------- */
  mat rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
  rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
  rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

  /* Replacing tt + 1 column on alpha_samples ------------- */
  if(augment_alpha){
    alpha_samples.col(tt + 1) =  alpha_samples.submat(index, tt_vec + 1);
  }

  if(partial){
    cube aug_rankings_index = aug_rankings.slices(index);
    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(span(0, num_obs - 1), span::all, span::all);
  }
}

void smc_mallows_new_users_reweight(
    vec& log_inc_wgt,
    rowvec& ESS_vec,
    vec& norm_wgt,
    const cube& aug_rankings,
    const mat& observed_rankings,
    const cube& rho_samples,
    const double& alpha,
    const mat& alpha_samples,
    const int& tt,
    const Rcpp::Nullable<vec> logz_estimate,
    const Rcpp::Nullable<vec> cardinalities,
    const int& num_obs,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& augment_alpha,
    const bool& partial,
    const std::string& metric = "footrule"
){
  int n_particles = rho_samples.n_rows;
  int n_items = rho_samples.n_cols;
  for (int ii{}; ii < n_particles; ++ii) {
    rowvec rho_samples_ii = \
      rho_samples(span(ii), span::all, span(tt + 1));

    double alpha_used = augment_alpha ? alpha_samples(ii, tt + 1) : alpha;


    /* Calculating log_z_alpha and log_likelihood ----------- */
    const double log_z_alpha = get_partition_function(
      n_items, alpha_used, cardinalities, logz_estimate, metric
    );

    double log_likelihood = get_exponent_sum(                          \
      alpha_used, rho_samples_ii.t(), n_items,
      partial ? aug_rankings(span(num_obs - num_new_obs, num_obs - 1), span::all, span(ii)) : observed_rankings,
      metric
    );
    log_inc_wgt(ii) = log_likelihood - num_new_obs * log_z_alpha - log(aug_prob(ii));
  }

  /* normalise weights ------------------------------------ */
  norm_wgt = normalize_weights(log_inc_wgt);

  ESS_vec(tt) = (sum(norm_wgt) * sum(norm_wgt)) / sum(norm_wgt % norm_wgt);
}
