#' @title Calculate Forward Probability
#' @description Function to calculate probability of assigning a set of specific ranks to an specific item
#' given its rank in the consensus ranking
#' @export
#'
#' @param item_ordering A vector of integer values to represent the specified queue of which unranked item to assign a rank for the proposed augmented ranking
#' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
#' @param remaining_set A vector of integer values to represent the elements (ranks) missing from original observed ranking
#' @param rho Numeric vector specifying the consensus ranking
#' @param alpha Numeric value og the scale parameter
#' @param n_items Integer is the number of items in a ranking
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @return List containing aug_ranking, a ranking sequence vector of the proposed augmented ranking and
#'         forward_prob a numerical value of the probability of creating the augmented ranking
#'         using the pseudolikelihood augmentation.
calculate_forward_probability = function(item_ordering, partial_ranking, remaining_set, rho, alpha, n_items, metric){


  # item ordering is the order of which items are assigned ranks in a specified order
  num_items_unranked = length(item_ordering)

  # prob of creating augmented ranking
  forward_auxiliary_ranking_probability = 1

  if(num_items_unranked == 1){

    # create new agumented ranking by sampling remaining ranks from set uniformly
    partial_ranking[is.na(partial_ranking)] <- remaining_set


  }else{

    auxiliary_ranking = rep(0, num_items_unranked)

    ########################################################
    ## LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY
    ########################################################
    # given the old and new item ordering and the list of missing rank, determine the sample probs for each iteration
    for (jj in 1:(num_items_unranked-1)){

      # items to sample rank
      item_to_sample_rank = item_ordering[jj]

      # the rank of item in rho
      rho_item_rank = rho[item_to_sample_rank]

      # next we get the sample probabilites for selecting a particular rank for an item
      # based on the current alpha and the rho rank for that item
      sample_prob_list = get_sample_probabilities(rho_item_rank = rho_item_rank, alpha = alpha,
                                                  remaining_set_ranks = remaining_set, metric = metric, n_items = n_items)
      #print(sample_prob_list)

      # fill in the new augmented ranking going forward
      auxiliary_ranking[jj] = sample(remaining_set, size = 1, replace = FALSE, prob = sample_prob_list)

      # save the probability of selecting the specific item rank in the old augmented ranking
      sample_prob = which(remaining_set == auxiliary_ranking[jj])
      forward_auxiliary_ranking_probability = forward_auxiliary_ranking_probability * (sample_prob_list[sample_prob])

      # remove selected auxiliary rank from the set of remaining possibles ranks to select
      remaining_set = remaining_set[ remaining_set != auxiliary_ranking[jj] ]

    }

    # last element in augmented ranking is deterministic - the prob is 1
    auxiliary_ranking[num_items_unranked] = remaining_set

    # fit the augmented ranking within the partial rankings with NAs
    partial_ranking[item_ordering] <- auxiliary_ranking # ranks for items

  } #end of if else statement


  output <- list("aug_ranking" = partial_ranking, "forward_prob" =  forward_auxiliary_ranking_probability)
  return(output)
}


#' @title Calculate Backward Probability
#' @description Function to calculate probability of assigning a set of specific ranks to an specific item
#' given its rank in the consensus ranking
#'
#' @param item_ordering A vector of integer values to represent the specified queue of which unranked item to assign a rank for the proposed augmented ranking
#' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
#' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obatined from calculate_forward_probability function
#' @param remaining_set A vector of integer values to represent the elements (ranks) missing from original observed ranking
#' @param rho Numeric vector specifying the consensus ranking
#' @param alpha Numeric value og the scale parameter
#' @param n_items Integer is the number of items in a ranking
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @return backward_auxiliary_ranking_probability A numerical value of creating the previous augmented ranking using the same item ordering used to create the
#' new auggmented ranking in calculate_forward_probability funtion.
#' @export
calculate_backward_probability = function(item_ordering, partial_ranking, current_ranking, remaining_set, rho, alpha, n_items, metric){
# given an old and new item ordering, sample ranking with new ordering and calc the forward and backward prob


  # item ordering is the order of which items are assigned ranks in a specified order
  num_items_unranked = length(item_ordering)
  num_ranks = length(current_ranking)

  # initialise prob of creating augmented ranking
  backward_auxiliary_ranking_probability = 1

  if(num_items_unranked != 1){

    # show the augmented parts of the current ranking
    current_ranking = current_ranking[item_ordering]

    ########################################################
    ## LOOP TO CALCULATE FORWARD AND BACKWARD PROBABILITY
    ########################################################
    # given the old and new item ordering and the list of missing rank, determine the sample probs for each iteration
    for (jj in 1:(num_items_unranked-1)){

      # items to sample rank
      item_to_sample_rank = item_ordering[jj]

      # the rank of item in rho
      rho_item_rank = rho[item_to_sample_rank]

      # next we get the sample probabilites for selecting a particular rank for an item
      # based on the current alpha and the rho rank for that item
      sample_prob_list = get_sample_probabilities(rho_item_rank = rho_item_rank, alpha = alpha,
                                                  remaining_set_ranks = remaining_set, metric = metric, n_items = n_items)

      # save the probability of selecting the specific item rank in the old augmented ranking
      sample_prob = which(remaining_set == current_ranking[jj])
      backward_auxiliary_ranking_probability = backward_auxiliary_ranking_probability * (sample_prob_list[sample_prob])

      # remove selected auxiliary rank from the set of remaining possibles ranks to select
      remaining_set = remaining_set[ remaining_set != current_ranking[jj] ]

    }
  }else{backward_auxiliary_ranking_probability = 1} # end of if else statement

  return(backward_auxiliary_ranking_probability)
}




#' @title Metropolis-Hastings Augmented Ranking (pseudolikelihood)
#' @description Function to perform Metropolis-Hastings for new augmented ranking using the pseudolikelihood augmentation approach
#'
#' @param alpha Numeric value og the scale parameter
#' @param rho Numeric vector specifying the consensus ranking
#' @param n_items Integer is the number of items in a ranking
#' @param partial_ranking An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
#' @param current_ranking An complete rank sequence vector of  the proposed augmented ranking obatined from calculate_forward_probability function
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @return = proposed augmented ranking or current ranking A ranking sequence vector representing proposed augmented ranking for next
#'         iteration of MCMC chain
#' @export
metropolis_hastings_aug_ranking_pseudo = function(alpha, rho, n_items, partial_ranking, current_ranking, metric){


  # augment incomplete ranks to initialise
  ranks = c(1:n_items)

  # find items missing from original observed ranking
  unranked_items = as.numeric(which(is.na(partial_ranking)))

  # find unallocated ranks from original observed ranking
  remaining_set = ranks[!ranks %in% partial_ranking]

  # if the observed and augmented ranking are exactly the same then break
  if(identical(partial_ranking,current_ranking) == TRUE){
    return(current_ranking)

  }else if(length(remaining_set)==1){

    return(current_ranking)

  }else{

    # randomly permute the unranked items to give the order in which they will be allocated
    item_ordering = sample(unranked_items, size=length(unranked_items), replace=F)
    proposal = calculate_forward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                             remaining_set = remaining_set, rho = rho, alpha = alpha,
                                             n_items = n_items, metric = metric)
    proposed_augmented_ranking = proposal$aug_ranking
    forward_prob = proposal$forward_prob

    backward_prob = calculate_backward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                                   current_ranking = current_ranking, remaining_set = remaining_set, rho = rho,
                                                   alpha = alpha, n_items = n_items, metric = metric)

    #############
    # MH TIME
    #############
    # Calculate the log posterior of the current and proposed rankings.
    # NB the current can usually be stored to save recalculating it, but we're not caring about that yet
    curr_logpost = get_mallows_loglik(alpha = alpha, rho = rho, n = n_items, rankings = t(current_ranking), metric = metric)
    prop_logpost = get_mallows_loglik(alpha = alpha, rho = rho, n = n_items, rankings = t(proposed_augmented_ranking), metric = metric)

    log_acceptance_prob = prop_logpost - curr_logpost - log(forward_prob) + log(backward_prob)

    if (log(runif(1))<log_acceptance_prob) {
      return(proposed_augmented_ranking)
    } else {
      return(current_ranking)
    }

  }
}




#' @title Correction Kernel (pseudolikelihood)
#' @description Function to determine if the augmented ranking is compatible with the new observed partial ranking.
#' If it is not, the we create a new augmentation using the pseudolikelihood approach and calculate the augmentation probability.
#'
#' @param R_obs An incomplete rank sequence vector of the original observed incomplete ranking which contains NAs
#' @param R_curr An complete rank sequence vector of  the proposed augmented ranking obatined from calculate_forward_probability function
#' @param rho Numeric vector specifying the consensus ranking
#' @param alpha Numeric value og the scale parameter
#' @param n_items Integer is the number of items in a ranking
#' @param metric A character string specifying the distance metric to use in the
#'   Bayesian Mallows Model. Available options are \code{"footrule"},
#'   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#'   \code{"ulam"}.
#' @return list containing R_obs, the proposed 'corrected' augmented ranking that is compatible with the new observed ranking for a user, and
#'         forward_auxiliary_ranking_probability, a numerical value for the probability of correcting the ranking to be compatible with R_obs.
#' @export
correction_kernel_pseudo = function(R_curr, R_obs, rho, alpha, n_items, metric){


  # check if new information means 'mistakes' made with augmented rankings
  check = (R_obs == R_curr)
  if (any(check==FALSE, na.rm = TRUE) == FALSE){

    forward_auxiliary_ranking_probability = 1
    output = list("ranking" = R_curr, "correction_prob" = forward_auxiliary_ranking_probability)

  }else{

    # resample from smaller pool of possible augmented rankings
    #select uniform the proposed ranking compatible with the known observed rankings
    ranks = c(1:n_items)

    # find items missing from original observed ranking
    unranked_items = as.numeric(which(is.na(R_obs)))

    # find elements missing from original observed ranking
    remaining_set = unique(ranks[!ranks %in% R_obs])

    # if we only have one missing rank, then we can
    if(length(unranked_items) == 1){

      # create new agumented ranking by sampling remaining ranks from set uniformly
      R_obs[is.na(R_obs)] <- remaining_set

      forward_auxiliary_ranking_probability = 1

      output = list("ranking" = R_obs, "correction_prob" = forward_auxiliary_ranking_probability )

    }else{

      # create new agumented ranking by using pseudo proposal
      item_ordering = sample(unranked_items, size=length(unranked_items), replace=F)
      # item ordering is the order of which items are assigned ranks in a specified order
      num_items_unranked = length(item_ordering)

      # creating now augmented ranking whilst simultaneously calculating the backwards prob of making the same
      # augmented ranking with an alternative item ordering
      auxiliary_ranking = rep(0, num_items_unranked)
      forward_auxiliary_ranking_probability = 1

      ########################################################
      ## Create new augmented ranking
      ########################################################
      # given the old and new item ordering and the list of missing rank, determine the sample probs for each iteration
      for (jj in 1:(num_items_unranked-1)){

        # items to sample rank
        item_to_sample_rank = item_ordering[jj]

        # the rank of item in rho
        rho_item_rank = rho[item_to_sample_rank]

        # next we get the sample probabilites for selecting a particular rank for an item
        # based on the current alpha and the rho rank for that item
        sample_prob_list = get_sample_probabilities(rho_item_rank = rho_item_rank, alpha = alpha,
                                                    remaining_set_ranks = remaining_set, metric = metric, n_items = n_items)

        # fill in the new augmented ranking going forward

        auxiliary_ranking[jj] = sample(remaining_set, size = 1, replace = FALSE, prob = sample_prob_list)

        sample_prob = which(remaining_set == auxiliary_ranking[jj])
        forward_auxiliary_ranking_probability = forward_auxiliary_ranking_probability * (sample_prob_list[sample_prob])

        # remove selected auxiliary rank from the set of remaining possibles ranks to select
        remaining_set = remaining_set[ remaining_set != auxiliary_ranking[jj] ]
      }

      # last element in augmented ranking is deterministic - the prob is 1
      auxiliary_ranking[num_items_unranked] = remaining_set

      # fit the augmented ranking within the partial rankings with NAs
      R_obs[item_ordering] <- auxiliary_ranking # ranks for items

      output = list("ranking" = R_obs, "correction_prob" = forward_auxiliary_ranking_probability )

    }

  }

  return(output)
}
