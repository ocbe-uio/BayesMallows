// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// asymptotic_partition_function
arma::vec asymptotic_partition_function(arma::vec alpha_vector, int n_items, std::string metric, int K, int n_iterations, double tol);
RcppExport SEXP _BayesMallows_asymptotic_partition_function(SEXP alpha_vectorSEXP, SEXP n_itemsSEXP, SEXP metricSEXP, SEXP KSEXP, SEXP n_iterationsSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha_vector(alpha_vectorSEXP);
    Rcpp::traits::input_parameter< int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type n_iterations(n_iterationsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(asymptotic_partition_function(alpha_vector, n_items, metric, K, n_iterations, tol));
    return rcpp_result_gen;
END_RCPP
}
// get_rank_distance
arma::vec get_rank_distance(arma::mat rankings, arma::vec rho, std::string metric);
RcppExport SEXP _BayesMallows_get_rank_distance(SEXP rankingsSEXP, SEXP rhoSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(get_rank_distance(rankings, rho, metric));
    return rcpp_result_gen;
END_RCPP
}
// compute_importance_sampling_estimate
arma::vec compute_importance_sampling_estimate(arma::vec alpha_vector, int n_items, std::string metric, int nmc);
RcppExport SEXP _BayesMallows_compute_importance_sampling_estimate(SEXP alpha_vectorSEXP, SEXP n_itemsSEXP, SEXP metricSEXP, SEXP nmcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type alpha_vector(alpha_vectorSEXP);
    Rcpp::traits::input_parameter< int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< int >::type nmc(nmcSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_importance_sampling_estimate(alpha_vector, n_items, metric, nmc));
    return rcpp_result_gen;
END_RCPP
}
// get_expected_distance
double get_expected_distance(double alpha, int n_items, std::string metric, const Rcpp::Nullable<arma::mat>& pfun_values);
RcppExport SEXP _BayesMallows_get_expected_distance(SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP metricSEXP, SEXP pfun_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::mat>& >::type pfun_values(pfun_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_expected_distance(alpha, n_items, metric, pfun_values));
    return rcpp_result_gen;
END_RCPP
}
// get_partition_function
double get_partition_function(double alpha, int n_items, std::string metric, const Rcpp::Nullable<arma::mat>& pfun_values);
RcppExport SEXP _BayesMallows_get_partition_function(SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP metricSEXP, SEXP pfun_valuesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::mat>& >::type pfun_values(pfun_valuesSEXP);
    rcpp_result_gen = Rcpp::wrap(get_partition_function(alpha, n_items, metric, pfun_values));
    return rcpp_result_gen;
END_RCPP
}
// rmallows
arma::mat rmallows(arma::vec rho0, double alpha0, int n_samples, int burnin, int thinning, int leap_size, std::string metric);
RcppExport SEXP _BayesMallows_rmallows(SEXP rho0SEXP, SEXP alpha0SEXP, SEXP n_samplesSEXP, SEXP burninSEXP, SEXP thinningSEXP, SEXP leap_sizeSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho0(rho0SEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< int >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(rmallows(rho0, alpha0, n_samples, burnin, thinning, leap_size, metric));
    return rcpp_result_gen;
END_RCPP
}
// run_mcmc
Rcpp::List run_mcmc(Rcpp::List data, Rcpp::List model_options, Rcpp::List compute_options, Rcpp::List priors, Rcpp::List initial_values, Rcpp::Nullable<arma::mat> pfun_values, Rcpp::Nullable<arma::mat> pfun_estimate, bool verbose);
RcppExport SEXP _BayesMallows_run_mcmc(SEXP dataSEXP, SEXP model_optionsSEXP, SEXP compute_optionsSEXP, SEXP priorsSEXP, SEXP initial_valuesSEXP, SEXP pfun_valuesSEXP, SEXP pfun_estimateSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type model_options(model_optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type compute_options(compute_optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initial_values(initial_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type pfun_values(pfun_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type pfun_estimate(pfun_estimateSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(run_mcmc(data, model_options, compute_options, priors, initial_values, pfun_values, pfun_estimate, verbose));
    return rcpp_result_gen;
END_RCPP
}
// run_smc
Rcpp::List run_smc(Rcpp::List data, Rcpp::List new_data, Rcpp::List model_options, Rcpp::List smc_options, Rcpp::List compute_options, Rcpp::List priors, Rcpp::List initial_values, Rcpp::Nullable<arma::mat> pfun_values, Rcpp::Nullable<arma::mat> pfun_estimate);
RcppExport SEXP _BayesMallows_run_smc(SEXP dataSEXP, SEXP new_dataSEXP, SEXP model_optionsSEXP, SEXP smc_optionsSEXP, SEXP compute_optionsSEXP, SEXP priorsSEXP, SEXP initial_valuesSEXP, SEXP pfun_valuesSEXP, SEXP pfun_estimateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type new_data(new_dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type model_options(model_optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type smc_options(smc_optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type compute_options(compute_optionsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type priors(priorsSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type initial_values(initial_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type pfun_values(pfun_valuesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type pfun_estimate(pfun_estimateSEXP);
    rcpp_result_gen = Rcpp::wrap(run_smc(data, new_data, model_options, smc_options, compute_options, priors, initial_values, pfun_values, pfun_estimate));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesMallows_asymptotic_partition_function", (DL_FUNC) &_BayesMallows_asymptotic_partition_function, 6},
    {"_BayesMallows_get_rank_distance", (DL_FUNC) &_BayesMallows_get_rank_distance, 3},
    {"_BayesMallows_compute_importance_sampling_estimate", (DL_FUNC) &_BayesMallows_compute_importance_sampling_estimate, 4},
    {"_BayesMallows_get_expected_distance", (DL_FUNC) &_BayesMallows_get_expected_distance, 4},
    {"_BayesMallows_get_partition_function", (DL_FUNC) &_BayesMallows_get_partition_function, 4},
    {"_BayesMallows_rmallows", (DL_FUNC) &_BayesMallows_rmallows, 7},
    {"_BayesMallows_run_mcmc", (DL_FUNC) &_BayesMallows_run_mcmc, 8},
    {"_BayesMallows_run_smc", (DL_FUNC) &_BayesMallows_run_smc, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesMallows(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
