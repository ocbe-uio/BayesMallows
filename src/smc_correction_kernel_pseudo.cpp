#include "RcppArmadillo.h"
#include "smc.h"
#include "misc.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
//' @title Correction Kernel (pseudolikelihood)
//' @description Function to determine if the augmented ranking is compatible with the new observed partial ranking.
//' If it is not, the we create a new augmentation using the pseudolikelihood approach and calculate the augmentation probability.
//'
//' @param observed_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
//' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obtained from calculate_forward_probability function
//' @param rho Numeric vector specifying the consensus ranking
//' @param alpha Numeric value of the scale parameter
//' @param n_items Integer is the number of items in a ranking
//' @param metric A character string specifying the distance metric to use in the
//'   Bayesian Mallows Model. Available options are \code{"footrule"},
//'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
//'   \code{"ulam"}.
//' @return list containing R_obs, the proposed 'corrected' augmented ranking that is compatible with the new observed ranking for a user, and
//'         forward_auxiliary_ranking_probability, a numerical value for the probability of correcting the ranking to be compatible with R_obs.
// [[Rcpp::export]]
Rcpp::List correction_kernel_pseudo(
    const arma::vec current_ranking,  //R_curr
    arma::vec observed_ranking, //R_obs
    const arma::vec rho,
    const double alpha,
    const int n_items,
    const std::string metric
) {
    bool observed_equals_current = arma::approx_equal(\
        observed_ranking, current_ranking, "absdiff", 0.1\
    );
    double correction_prob = 1.0;
    if (observed_equals_current) {
        return Rcpp::List::create(
            Rcpp::Named("ranking") = current_ranking,
            Rcpp::Named("correction_prob") = correction_prob
        );
    } else {
        // resample from smaller pool of possible augmented rankings
        // select uniform the proposed ranking compatible with the known observed rankings
        // ranks = c(1:n_items)

        // find items missing from original observed ranking
        const uvec& unranked_items = find_nonfinite(observed_ranking);

        // find elements missing from original observed ranking
        vec remaining_set = arma_setdiff_vec(current_ranking, observed_ranking);

        // if we only have one missing rank, then we can
        if (remaining_set.n_elem == 1) {

            // create new agumented ranking by sampling remaining ranks from set uniformly
            observed_ranking.elem(unranked_items) = remaining_set;
        } else {
            // create new agumented ranking by using pseudo proposal
            uvec item_ordering = new_pseudo_proposal(unranked_items);

            // item ordering is the order of which items are assigned ranks in a specified order
            const uword& num_items_unranked = item_ordering.n_elem;

            // creating now augmented ranking whilst simultaneously calculating the backwards prob of making the same
            // augmented ranking with an alternative item ordering
            ivec auxiliary_ranking = Rcpp::rep(0, num_items_unranked);

            //########################################################
            //## Create new augmented ranking
            //########################################################
            // given the item ordering and the list of missing rank, determine the sample probs for each iteration
            for (uword jj = 0; jj < num_items_unranked - 1; ++jj) {

                // items to sample rank
                const double& item_to_sample_rank = item_ordering(jj);

                // the rank of item in rho
                vec rho_item_rank;
                rho_item_rank = rho(item_to_sample_rank - 1);

                // next we get the sample probabilites for selecting a particular rank for
                // an item based on the current alpha and the rho rank for that item
                const vec sample_prob_list = get_sample_probabilities(\
                    rho_item_rank, alpha, remaining_set, metric, n_items\
                );

                // fill in the new augmented ranking going forward
                auxiliary_ranking(jj) = sample_one_with_prob(remaining_set, sample_prob_list);

                // save the probability of selecting the specific item rank in the old
                // augmented ranking
                const uvec& sample_prob = find(remaining_set == auxiliary_ranking(jj));
                correction_prob = \
                    correction_prob * \
                    arma::as_scalar(sample_prob_list(sample_prob));

                // remove selected auxiliary rank from the set of remaining possibles
                // ranks to select
                remaining_set = remaining_set(find(remaining_set != auxiliary_ranking(jj)));
            }
            // last element in augmented ranking is deterministic - the prob is 1
            auxiliary_ranking(num_items_unranked - 1) = arma::as_scalar(remaining_set);

            // fit the augmented ranking within the partial rankings with NAs
            vec ar;
            ar = conv_to<vec>::from(auxiliary_ranking);
            observed_ranking.elem(item_ordering - 1) = ar; // ranks for items
        }
        return Rcpp::List::create(
            Rcpp::Named("ranking") = observed_ranking,
            Rcpp::Named("correction_prob") = correction_prob
        );
    }
}
