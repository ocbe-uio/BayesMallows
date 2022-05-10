#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include <RcppArmadillo.h>

void smc_mallows_new_users_augment_partial(
    arma::cube& aug_rankings,
    arma::vec& aug_prob,
    const arma::cube rho_samples,
    const arma::mat& alpha_samples,
    const int& num_obs,
    const int& num_new_obs,
    const arma::mat& R_obs,
    const std::string& aug_method,
    const int& N,
    const std::string& metric,
    const int& tt,
    const int& n_items,
    const double& alpha,
    const bool& augment_alpha
);

void smc_mallows_new_users_reweight(
    arma::vec& log_inc_wgt,
    arma::rowvec& ESS_vec,
    arma::vec& norm_wgt,
    const arma::cube& aug_rankings,
    const arma::mat& observed_rankings,
    const arma::cube& rho_samples,
    const double& alpha,
    const arma::mat& alpha_samples,
    const int& N,
    const int& tt,
    const int& n_items,
    const Rcpp::Nullable<arma::vec> logz_estimate,
    const std::string& metric,
    const int& num_obs,
    const int& num_new_obs,
    const arma::vec& aug_prob,
    const bool& augment_alpha,
    const bool& partial
);

#endif
