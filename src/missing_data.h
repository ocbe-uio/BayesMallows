#ifndef MISSING_H
#define MISSING_H


void initialize_missing_ranks(arma::umat& rankings, const arma::umat& missing_indicator,
                              const arma::uvec& assessor_missing);

void update_missing_ranks(arma::umat& rankings, const arma::uvec& current_cluster_assignment,
                          arma::vec& aug_acceptance,
                          const arma::umat& missing_indicator,
                          const arma::uvec& assessor_missing,
                          const arma::vec& alpha, const arma::umat& rho,
                          const std::string& metric);

#endif
