#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include <RcppArmadillo.h>

void smc_mallows_new_users_augment_partial(arma::cube&, arma::vec&,
    const arma::mat&, const arma::vec&, const int&, const int&,
    const arma::mat&, const std::string&, const double&,
    const bool&, const std::string&);
void smc_mallows_new_users_reweight(
    arma::vec&, double&, arma::vec&,
    const arma::cube&, const arma::mat&, const arma::mat&, const double&,
    const arma::vec&, const Rcpp::Nullable<arma::vec>,
    const Rcpp::Nullable<arma::vec>,
    const int&, const int&, const arma::vec&, const bool&, const bool&,
    const std::string&);
void smc_mallows_new_users_resample(
    arma::mat&, arma::vec&, arma::cube&, const arma::vec&,
    const int& num_obs, const bool& augment_alpha, const bool& partial);
#endif
