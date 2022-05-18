// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// get_rank_distance
double get_rank_distance(arma::vec r1, arma::vec r2, std::string metric);
RcppExport SEXP _BayesMallows_get_rank_distance(SEXP r1SEXP, SEXP r2SEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type r1(r1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type r2(r2SEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(get_rank_distance(r1, r2, metric));
    return rcpp_result_gen;
END_RCPP
}
// rank_dist_sum
double rank_dist_sum(const arma::mat& rankings, const arma::vec& rho, const std::string& metric, const arma::vec& obs_freq);
RcppExport SEXP _BayesMallows_rank_dist_sum(SEXP rankingsSEXP, SEXP rhoSEXP, SEXP metricSEXP, SEXP obs_freqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type obs_freq(obs_freqSEXP);
    rcpp_result_gen = Rcpp::wrap(rank_dist_sum(rankings, rho, metric, obs_freq));
    return rcpp_result_gen;
END_RCPP
}
// rank_dist_vec
arma::vec rank_dist_vec(const arma::mat& rankings, const arma::vec& rho, const std::string& metric, const arma::vec& obs_freq);
RcppExport SEXP _BayesMallows_rank_dist_vec(SEXP rankingsSEXP, SEXP rhoSEXP, SEXP metricSEXP, SEXP obs_freqSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type obs_freq(obs_freqSEXP);
    rcpp_result_gen = Rcpp::wrap(rank_dist_vec(rankings, rho, metric, obs_freq));
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
// log_expected_dist
double log_expected_dist(const double& alpha, const int& n_items, const arma::vec& cardinalities, const std::string& metric);
RcppExport SEXP _BayesMallows_log_expected_dist(SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP cardinalitiesSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type cardinalities(cardinalitiesSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(log_expected_dist(alpha, n_items, cardinalities, metric));
    return rcpp_result_gen;
END_RCPP
}
// get_partition_function
double get_partition_function(int n_items, double alpha, const Rcpp::Nullable<arma::vec> cardinalities, const Rcpp::Nullable<arma::vec> logz_estimate, std::string metric);
RcppExport SEXP _BayesMallows_get_partition_function(SEXP n_itemsSEXP, SEXP alphaSEXP, SEXP cardinalitiesSEXP, SEXP logz_estimateSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type cardinalities(cardinalitiesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type logz_estimate(logz_estimateSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(get_partition_function(n_items, alpha, cardinalities, logz_estimate, metric));
    return rcpp_result_gen;
END_RCPP
}
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
// rmallows
arma::mat rmallows(arma::vec rho0, arma::vec obs_freq, double alpha0, int n_samples, int burnin, int thinning, int leap_size, std::string metric);
RcppExport SEXP _BayesMallows_rmallows(SEXP rho0SEXP, SEXP obs_freqSEXP, SEXP alpha0SEXP, SEXP n_samplesSEXP, SEXP burninSEXP, SEXP thinningSEXP, SEXP leap_sizeSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type rho0(rho0SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type obs_freq(obs_freqSEXP);
    Rcpp::traits::input_parameter< double >::type alpha0(alpha0SEXP);
    Rcpp::traits::input_parameter< int >::type n_samples(n_samplesSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type thinning(thinningSEXP);
    Rcpp::traits::input_parameter< int >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(rmallows(rho0, obs_freq, alpha0, n_samples, burnin, thinning, leap_size, metric));
    return rcpp_result_gen;
END_RCPP
}
// run_mcmc
Rcpp::List run_mcmc(arma::mat rankings, arma::vec obs_freq, int nmc, Rcpp::List constraints, Rcpp::Nullable<arma::vec> cardinalities, Rcpp::Nullable<arma::vec> logz_estimate, Rcpp::Nullable<arma::mat> rho_init, std::string metric, std::string error_model, int Lswap, int n_clusters, bool include_wcd, int leap_size, double alpha_prop_sd, double alpha_init, int alpha_jump, double lambda, double alpha_max, int psi, int rho_thinning, int aug_thinning, int clus_thin, bool save_aug, bool verbose, double kappa_1, double kappa_2, bool save_ind_clus);
RcppExport SEXP _BayesMallows_run_mcmc(SEXP rankingsSEXP, SEXP obs_freqSEXP, SEXP nmcSEXP, SEXP constraintsSEXP, SEXP cardinalitiesSEXP, SEXP logz_estimateSEXP, SEXP rho_initSEXP, SEXP metricSEXP, SEXP error_modelSEXP, SEXP LswapSEXP, SEXP n_clustersSEXP, SEXP include_wcdSEXP, SEXP leap_sizeSEXP, SEXP alpha_prop_sdSEXP, SEXP alpha_initSEXP, SEXP alpha_jumpSEXP, SEXP lambdaSEXP, SEXP alpha_maxSEXP, SEXP psiSEXP, SEXP rho_thinningSEXP, SEXP aug_thinningSEXP, SEXP clus_thinSEXP, SEXP save_augSEXP, SEXP verboseSEXP, SEXP kappa_1SEXP, SEXP kappa_2SEXP, SEXP save_ind_clusSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type obs_freq(obs_freqSEXP);
    Rcpp::traits::input_parameter< int >::type nmc(nmcSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type cardinalities(cardinalitiesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::vec> >::type logz_estimate(logz_estimateSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<arma::mat> >::type rho_init(rho_initSEXP);
    Rcpp::traits::input_parameter< std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< std::string >::type error_model(error_modelSEXP);
    Rcpp::traits::input_parameter< int >::type Lswap(LswapSEXP);
    Rcpp::traits::input_parameter< int >::type n_clusters(n_clustersSEXP);
    Rcpp::traits::input_parameter< bool >::type include_wcd(include_wcdSEXP);
    Rcpp::traits::input_parameter< int >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_init(alpha_initSEXP);
    Rcpp::traits::input_parameter< int >::type alpha_jump(alpha_jumpSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha_max(alpha_maxSEXP);
    Rcpp::traits::input_parameter< int >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< int >::type rho_thinning(rho_thinningSEXP);
    Rcpp::traits::input_parameter< int >::type aug_thinning(aug_thinningSEXP);
    Rcpp::traits::input_parameter< int >::type clus_thin(clus_thinSEXP);
    Rcpp::traits::input_parameter< bool >::type save_aug(save_augSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< double >::type kappa_1(kappa_1SEXP);
    Rcpp::traits::input_parameter< double >::type kappa_2(kappa_2SEXP);
    Rcpp::traits::input_parameter< bool >::type save_ind_clus(save_ind_clusSEXP);
    rcpp_result_gen = Rcpp::wrap(run_mcmc(rankings, obs_freq, nmc, constraints, cardinalities, logz_estimate, rho_init, metric, error_model, Lswap, n_clusters, include_wcd, leap_size, alpha_prop_sd, alpha_init, alpha_jump, lambda, alpha_max, psi, rho_thinning, aug_thinning, clus_thin, save_aug, verbose, kappa_1, kappa_2, save_ind_clus));
    return rcpp_result_gen;
END_RCPP
}
// calculate_backward_probability
double calculate_backward_probability(arma::uvec item_ordering, arma::vec partial_ranking, arma::vec current_ranking, arma::vec remaining_set, const arma::vec rho, const double alpha, const int n_items, const std::string metric);
RcppExport SEXP _BayesMallows_calculate_backward_probability(SEXP item_orderingSEXP, SEXP partial_rankingSEXP, SEXP current_rankingSEXP, SEXP remaining_setSEXP, SEXP rhoSEXP, SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type item_ordering(item_orderingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type partial_ranking(partial_rankingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type current_ranking(current_rankingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type remaining_set(remaining_setSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_backward_probability(item_ordering, partial_ranking, current_ranking, remaining_set, rho, alpha, n_items, metric));
    return rcpp_result_gen;
END_RCPP
}
// calculate_forward_probability
Rcpp::List calculate_forward_probability(arma::uvec item_ordering, arma::vec partial_ranking, arma::vec remaining_set, const arma::vec rho, const double alpha, const int n_items, const std::string metric);
RcppExport SEXP _BayesMallows_calculate_forward_probability(SEXP item_orderingSEXP, SEXP partial_rankingSEXP, SEXP remaining_setSEXP, SEXP rhoSEXP, SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type item_ordering(item_orderingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type partial_ranking(partial_rankingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type remaining_set(remaining_setSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(calculate_forward_probability(item_ordering, partial_ranking, remaining_set, rho, alpha, n_items, metric));
    return rcpp_result_gen;
END_RCPP
}
// correction_kernel
Rcpp::List correction_kernel(const arma::vec observed_ranking, const arma::vec current_ranking, const int n_items);
RcppExport SEXP _BayesMallows_correction_kernel(SEXP observed_rankingSEXP, SEXP current_rankingSEXP, SEXP n_itemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type observed_ranking(observed_rankingSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type current_ranking(current_rankingSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    rcpp_result_gen = Rcpp::wrap(correction_kernel(observed_ranking, current_ranking, n_items));
    return rcpp_result_gen;
END_RCPP
}
// correction_kernel_pseudo
Rcpp::List correction_kernel_pseudo(const arma::vec current_ranking, arma::vec observed_ranking, const arma::vec rho, const double alpha, const int n_items, const std::string metric);
RcppExport SEXP _BayesMallows_correction_kernel_pseudo(SEXP current_rankingSEXP, SEXP observed_rankingSEXP, SEXP rhoSEXP, SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type current_ranking(current_rankingSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type observed_ranking(observed_rankingSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(correction_kernel_pseudo(current_ranking, observed_ranking, rho, alpha, n_items, metric));
    return rcpp_result_gen;
END_RCPP
}
// get_exponent_sum
double get_exponent_sum(const double alpha, const arma::vec rho, const int n_items, arma::mat rankings, const std::string metric);
RcppExport SEXP _BayesMallows_get_exponent_sum(SEXP alphaSEXP, SEXP rhoSEXP, SEXP n_itemsSEXP, SEXP rankingsSEXP, SEXP metricSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    rcpp_result_gen = Rcpp::wrap(get_exponent_sum(alpha, rho, n_items, rankings, metric));
    return rcpp_result_gen;
END_RCPP
}
// get_sample_probabilities
arma::vec get_sample_probabilities(const arma::vec rho_item_rank, const double alpha, const arma::vec remaining_set_ranks, const std::string metric, const int n_items);
RcppExport SEXP _BayesMallows_get_sample_probabilities(SEXP rho_item_rankSEXP, SEXP alphaSEXP, SEXP remaining_set_ranksSEXP, SEXP metricSEXP, SEXP n_itemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type rho_item_rank(rho_item_rankSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type remaining_set_ranks(remaining_set_ranksSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    rcpp_result_gen = Rcpp::wrap(get_sample_probabilities(rho_item_rank, alpha, remaining_set_ranks, metric, n_items));
    return rcpp_result_gen;
END_RCPP
}
// leap_and_shift_probs
Rcpp::List leap_and_shift_probs(const arma::vec rho, const int leap_size, const int n_items);
RcppExport SEXP _BayesMallows_leap_and_shift_probs(SEXP rhoSEXP, SEXP leap_sizeSEXP, SEXP n_itemsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    rcpp_result_gen = Rcpp::wrap(leap_and_shift_probs(rho, leap_size, n_items));
    return rcpp_result_gen;
END_RCPP
}
// smc_mallows_new_item_rank
Rcpp::List smc_mallows_new_item_rank(const unsigned int& n_items, arma::cube& R_obs, const std::string& metric, const int& leap_size, const unsigned int& N, const unsigned int Time, const Rcpp::Nullable<arma::vec> logz_estimate, const int& mcmc_kernel_app, const double alpha, const double alpha_prop_sd, const double lambda, const double alpha_max, const std::string& aug_method, const bool verbose, const bool alpha_fixed);
RcppExport SEXP _BayesMallows_smc_mallows_new_item_rank(SEXP n_itemsSEXP, SEXP R_obsSEXP, SEXP metricSEXP, SEXP leap_sizeSEXP, SEXP NSEXP, SEXP TimeSEXP, SEXP logz_estimateSEXP, SEXP mcmc_kernel_appSEXP, SEXP alphaSEXP, SEXP alpha_prop_sdSEXP, SEXP lambdaSEXP, SEXP alpha_maxSEXP, SEXP aug_methodSEXP, SEXP verboseSEXP, SEXP alpha_fixedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const unsigned int& >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type R_obs(R_obsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const int& >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< const unsigned int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type logz_estimate(logz_estimateSEXP);
    Rcpp::traits::input_parameter< const int& >::type mcmc_kernel_app(mcmc_kernel_appSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_max(alpha_maxSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type aug_method(aug_methodSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type alpha_fixed(alpha_fixedSEXP);
    rcpp_result_gen = Rcpp::wrap(smc_mallows_new_item_rank(n_items, R_obs, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app, alpha, alpha_prop_sd, lambda, alpha_max, aug_method, verbose, alpha_fixed));
    return rcpp_result_gen;
END_RCPP
}
// smc_mallows_new_users
Rcpp::List smc_mallows_new_users(const arma::mat& R_obs, const std::string& type, const int& n_items, const std::string& metric, const int& leap_size, const int& N, int Time, const int& mcmc_kernel_app, const int& num_new_obs, const double alpha_prop_sd, const double lambda, const double alpha_max, const double alpha, const std::string& aug_method, const Rcpp::Nullable<arma::vec>& logz_estimate, const bool verbose);
RcppExport SEXP _BayesMallows_smc_mallows_new_users(SEXP R_obsSEXP, SEXP typeSEXP, SEXP n_itemsSEXP, SEXP metricSEXP, SEXP leap_sizeSEXP, SEXP NSEXP, SEXP TimeSEXP, SEXP mcmc_kernel_appSEXP, SEXP num_new_obsSEXP, SEXP alpha_prop_sdSEXP, SEXP lambdaSEXP, SEXP alpha_maxSEXP, SEXP alphaSEXP, SEXP aug_methodSEXP, SEXP logz_estimateSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type R_obs(R_obsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type type(typeSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const int& >::type leap_size(leap_sizeSEXP);
    Rcpp::traits::input_parameter< const int& >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type Time(TimeSEXP);
    Rcpp::traits::input_parameter< const int& >::type mcmc_kernel_app(mcmc_kernel_appSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_new_obs(num_new_obsSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_max(alpha_maxSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type aug_method(aug_methodSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec>& >::type logz_estimate(logz_estimateSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(smc_mallows_new_users(R_obs, type, n_items, metric, leap_size, N, Time, mcmc_kernel_app, num_new_obs, alpha_prop_sd, lambda, alpha_max, alpha, aug_method, logz_estimate, verbose));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_hastings_alpha
double metropolis_hastings_alpha(const double alpha, const int n_items, const arma::mat rankings, const std::string metric, const arma::vec rho, const Rcpp::Nullable<arma::vec> logz_estimate, const double alpha_prop_sd, const double lambda, const double alpha_max);
RcppExport SEXP _BayesMallows_metropolis_hastings_alpha(SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP rankingsSEXP, SEXP metricSEXP, SEXP rhoSEXP, SEXP logz_estimateSEXP, SEXP alpha_prop_sdSEXP, SEXP lambdaSEXP, SEXP alpha_maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const Rcpp::Nullable<arma::vec> >::type logz_estimate(logz_estimateSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_prop_sd(alpha_prop_sdSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha_max(alpha_maxSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_hastings_alpha(alpha, n_items, rankings, metric, rho, logz_estimate, alpha_prop_sd, lambda, alpha_max));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_hastings_aug_ranking
arma::vec metropolis_hastings_aug_ranking(const double& alpha, const arma::vec& rho, const int& n_items, const arma::vec& partial_ranking, const arma::vec& current_ranking, const std::string& metric, const bool& pseudo);
RcppExport SEXP _BayesMallows_metropolis_hastings_aug_ranking(SEXP alphaSEXP, SEXP rhoSEXP, SEXP n_itemsSEXP, SEXP partial_rankingSEXP, SEXP current_rankingSEXP, SEXP metricSEXP, SEXP pseudoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type partial_ranking(partial_rankingSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type current_ranking(current_rankingSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const bool& >::type pseudo(pseudoSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_hastings_aug_ranking(alpha, rho, n_items, partial_ranking, current_ranking, metric, pseudo));
    return rcpp_result_gen;
END_RCPP
}
// metropolis_hastings_rho
arma::vec metropolis_hastings_rho(const double alpha, const int n_items, const arma::mat rankings, const std::string metric, const arma::vec rho, const int leap_size);
RcppExport SEXP _BayesMallows_metropolis_hastings_rho(SEXP alphaSEXP, SEXP n_itemsSEXP, SEXP rankingsSEXP, SEXP metricSEXP, SEXP rhoSEXP, SEXP leap_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type n_items(n_itemsSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type rankings(rankingsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type metric(metricSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const int >::type leap_size(leap_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(metropolis_hastings_rho(alpha, n_items, rankings, metric, rho, leap_size));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_BayesMallows_get_rank_distance", (DL_FUNC) &_BayesMallows_get_rank_distance, 3},
    {"_BayesMallows_rank_dist_sum", (DL_FUNC) &_BayesMallows_rank_dist_sum, 4},
    {"_BayesMallows_rank_dist_vec", (DL_FUNC) &_BayesMallows_rank_dist_vec, 4},
    {"_BayesMallows_compute_importance_sampling_estimate", (DL_FUNC) &_BayesMallows_compute_importance_sampling_estimate, 4},
    {"_BayesMallows_log_expected_dist", (DL_FUNC) &_BayesMallows_log_expected_dist, 4},
    {"_BayesMallows_get_partition_function", (DL_FUNC) &_BayesMallows_get_partition_function, 5},
    {"_BayesMallows_asymptotic_partition_function", (DL_FUNC) &_BayesMallows_asymptotic_partition_function, 6},
    {"_BayesMallows_rmallows", (DL_FUNC) &_BayesMallows_rmallows, 8},
    {"_BayesMallows_run_mcmc", (DL_FUNC) &_BayesMallows_run_mcmc, 27},
    {"_BayesMallows_calculate_backward_probability", (DL_FUNC) &_BayesMallows_calculate_backward_probability, 8},
    {"_BayesMallows_calculate_forward_probability", (DL_FUNC) &_BayesMallows_calculate_forward_probability, 7},
    {"_BayesMallows_correction_kernel", (DL_FUNC) &_BayesMallows_correction_kernel, 3},
    {"_BayesMallows_correction_kernel_pseudo", (DL_FUNC) &_BayesMallows_correction_kernel_pseudo, 6},
    {"_BayesMallows_get_exponent_sum", (DL_FUNC) &_BayesMallows_get_exponent_sum, 5},
    {"_BayesMallows_get_sample_probabilities", (DL_FUNC) &_BayesMallows_get_sample_probabilities, 5},
    {"_BayesMallows_leap_and_shift_probs", (DL_FUNC) &_BayesMallows_leap_and_shift_probs, 3},
    {"_BayesMallows_smc_mallows_new_item_rank", (DL_FUNC) &_BayesMallows_smc_mallows_new_item_rank, 15},
    {"_BayesMallows_smc_mallows_new_users", (DL_FUNC) &_BayesMallows_smc_mallows_new_users, 16},
    {"_BayesMallows_metropolis_hastings_alpha", (DL_FUNC) &_BayesMallows_metropolis_hastings_alpha, 9},
    {"_BayesMallows_metropolis_hastings_aug_ranking", (DL_FUNC) &_BayesMallows_metropolis_hastings_aug_ranking, 7},
    {"_BayesMallows_metropolis_hastings_rho", (DL_FUNC) &_BayesMallows_metropolis_hastings_rho, 6},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesMallows(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
