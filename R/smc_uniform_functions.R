require("sets")
source("get_mallows_loglik.R")

metropolis_hastings_aug_ranking = function(R_curr, R_obs, alpha, rho, n_items, metric){

  # @description Function to perform Metropolis-Hastings for new augmented ranking

  # INPUT:
  #   @ param R_curr A ranking sequence vector of the current augmented ranking (no missing values)
  #   @param R_obs  A ranking sequence vector of the observed partial ranking (no missing values) The original incomplete partial ranking is in the rankings data set.
  #   @param alpha Numeric value og the scale parameter
  #   @param rho A ranking sequence vector
  #   @param n_items Integer is the number of items in a ranking
  #   @param metric A character string specifying the distance metric to use in the
  #   Bayesian Mallows Model. Available options are \code{"footrule"},
  #   \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
  #   \code{"ulam"}.

  # OUTPUT R_curr or R_obs A ranking sequence vector representing proposed augmented ranking for next
  #         iteration of MCMC chain

  #select uniform the proposed ranking compatible with the known observed rankings
  ranks = c(1:n_items)

  # find elements missing from original observed ranking
  remaining_set = unique(ranks[!ranks %in% R_obs])


  # if the observed and augmented ranking are exactly the same then break
  if(identical(R_obs,R_curr) == TRUE){
    print("identical")
    return(R_curr)

  }else if(length(remaining_set)==1){
    #print("one missing item rank only")
    return(R_curr)

  }else{

    # create new agumented ranking by sampling remaining ranks from set uniformly
    R_prop = R_obs
    R_prop[is.na(R_prop)] <- sample(remaining_set, size=length(remaining_set), replace=F)
    # if this doesn't work, make sure your ranking has NA in the vector, and not "NA"

    # evaluate the log-likelihood with current rankings
    mallows_lik_curr =  get_mallows_loglik(alpha = alpha, rho = rho, n = n_items, rankings = R_curr, metric = metric)
    mallows_lik_prop = get_mallows_loglik(alpha = alpha, rho = rho, n = n_items, rankings = R_prop, metric = metric)

    # calculate acceptance probability
    loga = mallows_lik_prop - mallows_lik_curr

    # determine whether to accept or reject proposed rho and return now consensus ranking
    p = runif(1, min = 0, max = 1)
    if(log(p) <= loga){
      return(R_prop)
    } else{
      return(R_curr)
    }

  }

}

correction_kernel <- function(R_curr, R_obs, n_items){

  # @description Function to determine if the augmented ranking is compatible with the new observed partial ranking.
  # If it is not, the we create a new augmentation using the random sampling approach and calculate the augmentation probability.

  # INPUT:
  #   @ param R_curr A ranking sequence vector of the current augmented ranking (no missing values)
  #   @param R_obs  A ranking sequence vector of the observed partial ranking (no missing values) The original incomplete partial ranking is in the rankings data set.
  #   @param n_items Integer is the number of items in a ranking

  # OUTPUT: List containing the proposed 'corrected' augmented ranking that is compatible with the new observed ranking for a user


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
    if(length(remaining_set == 1)){R_prop[is.na(R_prop)] = remaining_set}
    else {R_prop[is.na(R_prop)] <- sample(remaining_set, size=sum(is.na(R_prop)), replace=F)}

    correction_prob = 1/factorial(length(remaining_set))
    output = list("ranking" = R_prop, "correction_prob" = correction_prob )
    return(output)
  }

}



###############
# test script
###############
# # start a new session of R before running this script!
#
# # This functions uses get_mallows_log_lik so if the checks match in the worker function
# # then it is very likely that this function will return the correct outputs.
#


set.seed(101)
require("BayesMallows")

###############################################
# tests for M-H_aug_ranking function
###############################################


rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6


R_curr = c(1,2,3,6,5,4)
R_obs = c(1,2,3,NA,NA,NA)
test_1 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric)
print(test_1)
# 1 2 3 4 6 5
# if rho != rho_prime, then it should have a ulam distance of 1
# if rho == rho_prime, then it should have ulam distance of 0
BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
# should be 1


# R_curr == rho so we should get rho
R_curr = rho
R_obs = c(1,2,3,NA,NA,NA)
test_2 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric)
print(test_2)
# 1 2 3 4 5 6

all(test_2 == rho)
# should be TRUE

# if rho != rho_prime, then it should have a ulam distance of 1
# if rho == rho_prime, then it should have ulam distance of 0
BayesMallows:::get_rank_distance(rho, test_2, metric= "ulam")
# should be 0


# R_curr == rho so we should get rho

R_curr = c(1,2,3,6,5,4)
R_obs = c(1,2,3,6,5,NA)
test_3 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric)
print(test_3)
# 1 2 3 6 5 4

all(test_3 == R_curr)
# should be TRUE




#########################################################
### tests relating to the correction_kernel function
#########################################################


set.seed(101)
n_items= 6

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,NA,NA,NA)

test_4 = correction_kernel(R_curr, R_obs, n_items)
print(test_4)
#$ranking
#[1] 1 2 3 4 5 6
#
#$correction_prob
#[1] 1


all(test_4$ranking == R_curr) == TRUE
# TRUE


R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,5,NA,NA)
test_5 = correction_kernel(R_curr, R_obs, n_items)
print(test_5)
#$ranking
#[1] 1 2 3 5 4 6
#
#$correction_prob
#[1] 0.5


all(test_5$ranking == R_curr) == TRUE
#  should be FALSE


R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,4,5,6)
test_6 = correction_kernel(R_curr, R_obs, n_items)
print(test_6)
#$ranking
#[1] 1 2 3 4 5 6
#
#$correction_prob
#[1] 1


all(test_6$ranking == R_curr) == TRUE
#  should be TRUE










