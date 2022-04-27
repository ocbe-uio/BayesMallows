#include <RcppArmadillo.h>
#include "leapandshift.h"
#include "distances.h"

// [[Rcpp::depends(RcppArmadillo)]]

// Update shape parameters for the Bernoulli error model
void update_shape_bernoulli(
    double& shape_1,
    double& shape_2,
    const double& kappa_1,
    const double& kappa_2,
    const arma::mat& rankings,
    const Rcpp::List& constraints
){

  int n_items = rankings.n_rows;
  int n_assessors = rankings.n_cols;
  int sum_1 = 0, sum_2 = 0;
  for(int i = 0; i < n_assessors; ++i){

    Rcpp::List assessor_constraints = Rcpp::as<Rcpp::List>(constraints[i]);
    for(int j = 0; j < n_items; ++j){

      arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[j]);
      arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[j]);

      for(unsigned int k = 0; k < items_above.size(); ++k){
        int g = (arma::as_scalar(rankings.col(i).row(j)) < arma::as_scalar(rankings.col(i).row(items_above(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
      for(unsigned int k = 0; k < items_below.size(); ++k){
        int g = (arma::as_scalar(rankings.col(i).row(j)) > arma::as_scalar(rankings.col(i).row(items_below(k) - 1)));
        sum_1 += g;
        sum_2 += 1 - g;
      }
    }
  }
  shape_1 = kappa_1 + sum_1;
  shape_2 = kappa_2 + sum_2;

}

void find_pairwise_limits(int& left_limit, int& right_limit, const int& item,
                          const Rcpp::List& assessor_constraints,
                          const arma::vec& current_ranking){

  // Find the items which are preferred to the given item
  // Items go from 1, ..., n_items, so must use [item - 1]
  arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[item - 1]);
  arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[item - 1]);

  // If there are any items above, we must find the possible rankings
  if(items_above.n_elem > 0){
    // Again subtracting 1 because of zero-first indexing
    // Find all the rankings of the items that are preferred to *item*
    arma::vec rankings_above = current_ranking.elem(items_above - 1);
    left_limit = arma::max(rankings_above);
  }

  // If there are any items below, we must find the possible rankings
  if(items_below.n_elem > 0){
    // Find all the rankings of the items that are disfavored to *item*
    arma::vec rankings_below = current_ranking.elem(items_below - 1);
    right_limit = arma::min(rankings_below);
  }

}

arma::vec propose_pairwise_augmentation(const arma::vec& ranking, const Rcpp::List& assessor_constraints){

  int n_items = ranking.n_elem;

  // Extract the constraints for this particular assessor
  arma::uvec constrained_items = Rcpp::as<arma::uvec>(assessor_constraints[0]);

  // Sample an integer between 1 and n_items
  int item = arma::randi<int>(arma::distr_param(1, n_items));

  // Left and right limits of the interval we draw ranks from
  // Correspond to l_j and r_j, respectively, in Vitelli et al. (2018), JMLR, Sec. 4.2.
  int left_limit = 0, right_limit = n_items + 1;
  find_pairwise_limits(left_limit, right_limit, item, assessor_constraints, ranking);

  // Now complete the leap step by sampling a new proposal uniformly between
  // left_limit + 1 and right_limit - 1
  int proposed_rank = arma::randi<int>(arma::distr_param(left_limit + 1, right_limit - 1));

  // Assign the proposal to the (item-1)th item
  arma::vec proposal = ranking;
  proposal(item - 1) = proposed_rank;

  double delta_r;
  arma::uvec indices;

  // Do the shift step
  shift_step(proposal, ranking, item, delta_r, indices);

  return proposal;
}

arma::vec propose_swap(const arma::vec& ranking, const Rcpp::List& assessor_constraints,
                       int& g_diff, const int& Lswap){

  int n_items = ranking.n_elem;

  // Draw a random number, representing an item
  int u = arma::as_scalar(arma::randi(1, arma::distr_param(1, n_items - Lswap)));

  int ind1 = arma::as_scalar(arma::find(ranking == u));
  int ind2 = arma::as_scalar(arma::find(ranking == (u + Lswap)));

  arma::vec proposal = ranking;
  proposal(ind1) = ranking(ind2);
  proposal(ind2) = ranking(ind1);

  // First consider the first item that was switched
  arma::uvec items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[ind1]);
  arma::uvec items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[ind1]);

  for(unsigned int j = 0; j < items_above.size(); ++j){
    g_diff += (proposal(items_above[j] - 1) > proposal(ind1)) - (ranking(items_above[j] - 1) > ranking(ind1));
  }
  for(unsigned int j = 0; j < items_below.size(); ++j){
    g_diff += (proposal(items_below[j] - 1) < proposal(ind1)) - (ranking(items_below[j] - 1) < ranking(ind1));
  }

  // Now consider the second item
  items_above = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[1])[ind2]);
  items_below = Rcpp::as<arma::uvec>(Rcpp::as<Rcpp::List>(assessor_constraints[2])[ind2]);

  for(unsigned int j = 0; j < items_above.size(); ++j){
    g_diff += (proposal(items_above[j] - 1) > proposal(ind1)) - (ranking(items_above[j] - 1) > ranking(ind1));
  }
  for(unsigned int j = 0; j < items_below.size(); ++j){
    g_diff += (proposal(items_below[j] - 1) < proposal(ind1)) - (ranking(items_below[j] - 1) < ranking(ind1));
  }
  return proposal;
}


void augment_pairwise(
    arma::mat& rankings,
    const arma::uvec& current_cluster_assignment,
    const arma::vec& alpha,
    const double& theta,
    const arma::mat& rho,
    const std::string& metric,
    const Rcpp::List& constraints,
    arma::vec& aug_acceptance,
    const bool& clustering,
    const std::string& error_model,
    const int& Lswap
){

  int n_assessors = rankings.n_cols;
  int n_items = rankings.n_rows;

  for(int i = 0; i < n_assessors; ++i){

    arma::vec proposal;
    // Summed difference over error function before and after proposal
    int g_diff = 0;

    // Sample a proposal, depending on the error model
    if(error_model == "none"){
      proposal = propose_pairwise_augmentation(rankings.col(i), Rcpp::as<Rcpp::List>(constraints[i]));
    } else if(error_model == "bernoulli"){
      proposal = propose_swap(rankings.col(i), Rcpp::as<Rcpp::List>(constraints[i]), g_diff, Lswap);
    } else {
      Rcpp::stop("error_model must be 'none' or 'bernoulli'");
    }

    // Finally, decide whether to accept the proposal or not
    // Draw a uniform random number
    double u = std::log(arma::randu<double>());

    // Find which cluster the assessor belongs to
    int cluster = current_cluster_assignment(i);

    double ratio = -alpha(cluster) / n_items *
      (get_rank_distance(proposal, rho.col(cluster), metric) -
      get_rank_distance(rankings.col(i), rho.col(cluster), metric));

    if((theta > 0) & (g_diff != 0)) {
      ratio += g_diff * std::log(theta / (1 - theta));
    }

    if(ratio > u){
      rankings.col(i) = proposal;
      ++aug_acceptance(i);
    }
  }

}

