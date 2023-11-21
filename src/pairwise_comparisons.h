#ifndef PAIRWISE_H
#define PAIRWISE_H


void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const double& theta,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    const std::string& error_model,
    const unsigned int& swap_leap
);


#endif
