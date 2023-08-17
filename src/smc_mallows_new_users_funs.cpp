#include <RcppArmadillo.h>
#include "misc.h"
#include "sample.h"
#include "setdiff.h"
#include "smc.h"
#include "partitionfuns.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

void smc_mallows_new_users_augment_partial(
    arma::ucube& aug_rankings,
    arma::vec& aug_prob,
    const arma::ucube& rho_samples,
    const arma::mat& alpha_samples,
    const uint& num_obs,
    const uint& num_new_obs,
    const arma::umat& R_obs,
    const std::string& aug_method,
    const int& tt,
    const double& alpha,
    const bool& augment_alpha,
    const std::string& metric = "footrule"
){
  uint N = rho_samples.n_rows;
  uint n_items = rho_samples.n_cols;
  uvec ranks = regspace<uvec>(1, n_items);
  for (uint ii{}; ii < N; ++ii) {
    for (uint jj = num_obs - num_new_obs; jj < num_obs; ++jj) {
      uvec partial_ranking = R_obs.row(jj).t();

      // find items missing from original observed ranking
      const uvec& unranked_items = find(partial_ranking == 0);

      // find ranks missing from ranking
      const uvec& missing_ranks = conv_to<uvec>::from(setdiff_template(ranks, partial_ranking));

      // fill in missing ranks based on choice of augmentation method
      Rcpp::List proposal{};
      const bool pseudo = is_pseudo(aug_method, metric);
      if (!pseudo) {
        // create new augmented ranking by sampling remaining ranks from set uniformly
        partial_ranking.elem(find(partial_ranking == 0)) = shuffle(missing_ranks);

        aug_rankings(span(jj), span::all, span(ii)) = partial_ranking;
        aug_prob(ii) = divide_by_fact(aug_prob(ii), missing_ranks.n_elem);
      } else {
        // randomly permute the unranked items to give the order in which they will be allocated
        uvec item_ordering = shuffle(unranked_items);
        const urowvec& rho_s = rho_samples(span(ii), span::all, span(tt + 1));
        proposal = calculate_forward_probability(
          item_ordering, partial_ranking, missing_ranks, rho_s.t(),
          augment_alpha ? alpha_samples(ii, tt + 1) : alpha, n_items, metric
        );

        const uvec& a_rank = proposal["aug_ranking"];
        const double& f_prob = proposal["forward_prob"];
        aug_rankings(span(jj), span::all, span(ii)) = a_rank;
        aug_prob(ii) = aug_prob(ii) * f_prob;
      }
    }
  }
}

vec initialize_alpha(const int& N){
  return Rcpp::rexp(N, 1);
}

void smc_mallows_new_users_resample(
    ucube& rho_samples, mat& alpha_samples, ucube& aug_rankings,
    const vec& norm_wgt, const int& tt, const int& num_obs,
    const bool& augment_alpha, const bool& partial
){
  int N = rho_samples.n_rows;
  /* Resample particles using multinomial resampling ------ */
  uvec index = sample(regspace<uvec>(0, N - 1), N, true, norm_wgt);
  uvec tt_vec;
  tt_vec = tt;

  /* Replacing tt + 1 slice on rho_samples ---------------- */
  umat rho_samples_slice_11p1 = rho_samples.slice(tt + 1);
  rho_samples_slice_11p1 = rho_samples_slice_11p1.rows(index);
  rho_samples.slice(tt + 1) = rho_samples_slice_11p1;

  /* Replacing tt + 1 column on alpha_samples ------------- */
  if(augment_alpha){
    alpha_samples.col(tt + 1) =  alpha_samples.submat(index, tt_vec + 1);
  }

  if(partial){
    ucube aug_rankings_index = aug_rankings.slices(index);
    aug_rankings.rows(0, num_obs - 1) = aug_rankings_index(span(0, num_obs - 1), span::all, span::all);
  }
}

void smc_mallows_new_users_reweight(
    vec& log_inc_wgt,
    rowvec& ESS_vec,
    vec& norm_wgt,
    const ucube& aug_rankings,
    const umat& observed_rankings,
    const ucube& rho_samples,
    const double& alpha,
    const mat& alpha_samples,
    const int& tt,
    const Rcpp::Nullable<vec> logz_estimate,
    const int& num_obs,
    const int& num_new_obs,
    const vec& aug_prob,
    const bool& augment_alpha,
    const bool& partial,
    const std::string& metric = "footrule"
){
  int N = rho_samples.n_rows;
  uint n_items = rho_samples.n_cols;
  for (int ii{}; ii < N; ++ii) {
    Rcpp::Nullable<vec> cardinalities = R_NilValue;

    urowvec rho_samples_ii = \
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
