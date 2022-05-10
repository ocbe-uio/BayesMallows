#ifndef NEW_USERS_AUG
#define NEW_USERS_AUG

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
#endif
