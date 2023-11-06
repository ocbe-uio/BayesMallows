#ifndef MISSING_H
#define MISSING_H


void initialize_missing_ranks(arma::mat& rankings, const arma::umat& missing_indicator);

void update_missing_ranks(arma::mat& rankings, const arma::uvec& current_cluster_assignment,
                          const arma::umat& missing_indicator,
                          const arma::vec& alpha, const arma::mat& rho,
                          const std::string& metric);

#endif
