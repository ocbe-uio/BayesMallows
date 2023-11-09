#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include "RcppArmadillo.h"

void smc_mallows_new_users_augment_partial(arma::cube&, arma::vec&,
    const arma::mat&, const arma::vec&, const int&,
    const std::string&,
    const arma::umat& missing_indicator, const std::string&);
void smc_mallows_new_users_reweight(
    double&, arma::vec&,
    const arma::cube&, const arma::mat&, const arma::mat&,
    const arma::vec&, const Rcpp::List&,
    const int&, const arma::vec&, const bool&,
    const std::string&);
void smc_mallows_new_users_resample(
    arma::mat&, arma::vec&, arma::cube&, const arma::vec&,
    const bool& partial);


#endif
