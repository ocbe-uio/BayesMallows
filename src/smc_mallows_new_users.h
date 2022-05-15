#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include <RcppArmadillo.h>

arma::mat initialize_alpha(const int& N, const int& Time);
void smc_mallows_new_users_augment_partial(arma::cube&, arma::vec&,
    const arma::cube&, const arma::mat&, const int&, const int&,
    const arma::mat&, const std::string&, const std::string&, const int&,
    const double&, const bool&);
void smc_mallows_new_users_reweight(
    arma::vec&, arma::rowvec&, arma::vec&,
    const arma::cube&, const arma::mat&, const arma::cube&, const double&,
    const arma::mat&, const int&, const Rcpp::Nullable<arma::vec>, const std::string&,
    const int&, const int&, const arma::vec&, const bool&, const bool&);
void smc_mallows_new_users_resample(
    arma::cube&, arma::mat&, arma::cube&, const arma::vec&, const int& tt,
    const int& num_obs, const bool& augment_alpha, const bool& partial);
#endif
