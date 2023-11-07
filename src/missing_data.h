#ifndef MISSING_H
#define MISSING_H

arma::vec propose_augmentation(const arma::vec& ranks, const arma::uvec& indicator);

arma::vec make_new_augmentation(const arma::vec& rankings, const arma::uvec& missing_indicator,
                                const double& alpha, const arma::vec& rho,
                                const std::string& metric);

void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator);

void update_missing_ranks(arma::mat& rankings, const arma::uvec& current_cluster_assignment,
                          const arma::umat& missing_indicator,
                          const arma::vec& alpha, const arma::mat& rho,
                          const std::string& metric);

#endif
