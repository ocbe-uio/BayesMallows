#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include <RcppArmadillo.h>

void smc_mallows_new_users_augment_partial(arma::ucube&, arma::vec&,
    const arma::ucube&, const arma::mat&, const uint&, const uint&,
    const arma::umat&, const std::string&, const int&, const double&,
    const bool&, const std::string&);
void smc_mallows_new_users_reweight(
    arma::vec&, arma::rowvec&, arma::vec&,
    const arma::ucube&, const arma::umat&, const arma::ucube&, const double&,
    const arma::mat&, const int&, const Rcpp::Nullable<arma::vec>,
    const int&, const int&, const arma::vec&, const bool&, const bool&,
    const std::string&);
void smc_mallows_new_users_resample(
    arma::ucube&, arma::mat&, arma::ucube&, const arma::vec&, const int& tt,
    const int& num_obs, const bool& augment_alpha, const bool& partial);
#endif
