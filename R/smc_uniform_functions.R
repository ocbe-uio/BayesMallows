#' @title Correction Kernel
#' @description Function to determine if the augmented ranking is compatible with the new observed partial ranking.
#' If it is not, the we create a new augmentation using the random sampling approach and calculate the augmentation probability.
#'
#' @param R_curr A ranking sequence vector of the current augmented ranking (no missing values)
#' @param R_obs  A ranking sequence vector of the observed partial ranking (no missing values) The original incomplete partial ranking is in the rankings data set.
#' @param n_items Integer is the number of items in a ranking
#' @export
#' @return List containing the proposed 'corrected' augmented ranking that is compatible with the new observed ranking for a user
correction_kernel <- function(R_curr, R_obs, n_items){



  # check if new information means 'mistakes' made with augmented rankings
  check = (R_obs == R_curr)
  if (any(check==FALSE, na.rm = TRUE) == FALSE){
    correction_prob = 1
    output = list("ranking" = R_curr, "correction_prob" = correction_prob)
    return(output)

  } else {

    # resample from smaller pool of possible augmented rankings
    #select uniform the proposed ranking compatible with the known observed rankings
    ranks = c(1:n_items)

    # find elements missing from original observed ranking
    remaining_set = unique(ranks[!ranks %in% R_obs])

    # create new agumented ranking by sampling remaining ranks from set uniformly
    R_prop = R_obs
    if(length(remaining_set == 1)){
      R_prop[is.na(R_prop)] = remaining_set
    } else {
      R_prop[is.na(R_prop)] <- sample(remaining_set, size=sum(is.na(R_prop)), replace=F)
    }

    correction_prob = 1/factorial(length(remaining_set))
    output = list("ranking" = R_prop, "correction_prob" = correction_prob )
    return(output)
  }

}
