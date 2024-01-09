#include "classes.h"
#include "leapandshift.h"
#include "proposal_functions.h"
using namespace arma;

AlphaRatio make_new_alpha(
    const double& alpha_old,
    const vec& rho_old,
    const double& alpha_prop_sd,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartitionFunction>& pfun,
    const arma::mat& rankings,
    const arma::vec& observation_frequency,
    const double& n_items,
    const Priors& priors) {

  double alpha_proposal = R::rlnorm(std::log(alpha_old), alpha_prop_sd);
  double rank_dist = sum(distfun->d(rankings, rho_old) % observation_frequency);
  double alpha_diff = alpha_old - alpha_proposal;

  double ratio =
    alpha_diff / n_items * rank_dist +
    priors.lambda * alpha_diff +
    sum(observation_frequency) *
    (pfun->logz(alpha_old) - pfun->logz(alpha_proposal)) +
    std::log(alpha_proposal) - std::log(alpha_old);

  return AlphaRatio{alpha_proposal, ratio > std::log(R::unif_rand())};
}

vec make_new_rho(
    const vec& current_rho,
    const mat& rankings,
    double alpha_old,
    int leap_size,
    const std::unique_ptr<Distance>& distfun,
    vec observation_frequency) {

  vec rho_proposal;
  uvec indices;
  double prob_backward, prob_forward;
  int n_items = current_rho.n_elem;

  leap_and_shift(
    rho_proposal, indices, prob_backward, prob_forward,
    current_rho, leap_size, distfun);

  const arma::mat& r = rankings.rows(indices);
  double dist_new = arma::sum(
    distfun->d(r, rho_proposal(indices)) % observation_frequency
  );
  double dist_old = arma::sum(
    distfun->d(r, current_rho(indices)) % observation_frequency
  );

  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(prob_backward) - std::log(prob_forward);

  if(ratio > std::log(R::unif_rand())){
    return(rho_proposal);
  } else {
    return(current_rho);
  }
}



