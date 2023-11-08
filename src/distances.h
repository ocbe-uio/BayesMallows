#ifndef DISTANCES_H
#define DISTANCES_H



double get_rank_distance(arma::vec, arma::vec, std::string);
double rank_dist_sum(const arma::mat&, const arma::vec&, const std::string&, const arma::vec&);

arma::vec rank_dist_vec(const arma::mat& rankings, const arma::vec& rho,
                        const std::string& metric, const arma::vec& observation_frequency);

#endif
