#include "classes.h"
#include "leapandshift.h"
#include "proposal_functions.h"
using namespace arma;

AlphaRatio make_new_alpha(
    double alpha_old,
    const vec& rho_old,
    double alpha_prop_sd,
    const std::unique_ptr<Distance>& distfun,
    const std::unique_ptr<PartitionFunction>& pfun,
    const arma::mat& rankings,
    const arma::vec& observation_frequency,
    double n_items,
    const Priors& priors) {

  double alpha_proposal = R::rlnorm(std::log(alpha_old), alpha_prop_sd);
  double rank_dist = sum(distfun->matdist(rankings, rho_old) % observation_frequency);
  double alpha_diff = alpha_old - alpha_proposal;

  double ratio =
    alpha_diff / n_items * rank_dist +
    priors.lambda * alpha_diff +
    sum(observation_frequency) *
    (pfun->logz(alpha_old) - pfun->logz(alpha_proposal)) +
    priors.gamma * (std::log(alpha_proposal) - std::log(alpha_old));

  return AlphaRatio{alpha_proposal, ratio > std::log(R::unif_rand())};
}

vec make_new_rho(
    const vec& current_rho,
    const mat& rankings,
    double alpha_old,
    int leap_size,
    const std::unique_ptr<Distance>& distfun,
    vec observation_frequency) {

  int n_items = current_rho.n_elem;

  LeapShiftObject ls = leap_and_shift(current_rho, leap_size, distfun);

  double dist_new = arma::sum(
    distfun->matdist(rankings.rows(ls.indices), ls.rho_proposal(ls.indices)) % observation_frequency
  );
  double dist_old = arma::sum(
    distfun->matdist(rankings.rows(ls.indices), current_rho(ls.indices)) % observation_frequency
  );

  double ratio = - alpha_old / n_items * (dist_new - dist_old) +
    std::log(ls.prob_backward) - std::log(ls.prob_forward);

  if(ratio > std::log(R::unif_rand())){
    return(ls.rho_proposal);
  } else {
    return(current_rho);
  }
}



