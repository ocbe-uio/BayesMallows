#ifndef SMC_NEW_USERS
#define SMC_NEW_USERS

#include "RcppArmadillo.h"
#include "classes.h"

void reweight_new_users(
    SMCParameters& pars,
    const SMCAugmentation& aug,
    const SMCData& dat,
    const Rcpp::List& logz_list
);


#endif
