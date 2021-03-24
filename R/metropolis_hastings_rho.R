metropolis_hastings_rho <- function(alpha, n_items, rankings, metric, rho, leap_size){

  # @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!

  # INPUT:
  #   @param alpha Numeric value og the scale parameter
  #   @param n_items Integer is the number of items in a ranking
  #   @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
  #   rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
  #   can be a vector.
  #   @param rho A ranking sequence vector
  #   @param leap_size Integer specifying the step size of the leap-and-shift
  #   proposal distribution.

  # OUTPUT: rho or rho_prime A ranking sequence vector to be the next value of rho in the MCMC chain

  # create new potential consensus ranking
  kernel <- leap_and_shift_probs(rho=rho , n_items=n_items, leap_size=leap_size)
  print(kernel)

  # output from leap-and-shift is of the following
  #leap_shift_list <- list("rho_prime" = rho_prime, "forwards_prob" = forwards_prob, "backwards_prob" = backwards_prob)
  rho_prime = kernel$rho_prime
  forwards_prob = kernel$forwards_prob #rho_prime|rho
  backwards_prob = kernel$backwards_prob #rho|rho_prime

  # evaluate the log-likelihood with current rankings
  mallows_loglik_curr =  get_mallows_loglik(alpha = alpha, rho = rho, n_items = n_items, rankings = rankings, metric = metric)
  print(mallows_loglik_curr)
  mallows_loglik_prop = get_mallows_loglik(alpha = alpha, rho = rho_prime, n_items = n_items, rankings = rankings, metric = metric)
  print(mallows_loglik_prop)

  # calculate acceptance probability
  loga =  log(backwards_prob) - log(forwards_prob) + mallows_loglik_prop - mallows_loglik_curr


  # determine whether to accept or reject proposed rho and return now consensus ranking
  p = runif(1, min = 0, max = 1)
  if(log(p) <= loga){
    return(rho_prime)
  } else{
    return(rho)
  }

}

###############
# test script
###############
# # start a new session of R before running this script!
#
# # This functions uses get_mallows_log_lik and leap_and_shift_probs so if the checks match in those worker functions
# # then it is very likely that this function will return the correct outputs.
#
# set.seed(101)
#
# rho = c(1,2,3,4,5,6)
# alpha = 2
# metric = "footrule"
# n_items= 6
#
# rankings =   sample_mallows(rho0 = rho, alpha0 = alpha, n_samples = 10,
#                             burnin = 1000, thinning = 500)
#
# # you can confirm the print statements inside the metropolis_hastings_rho match get_mallows_loglik and leap_and_shift_probs
# test_1 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 1)
# print(test_1)
# # 1 2 3 5 4 6
# # if rho != rho_prime, then it should have a ulam distance of 1
# # if rho == rho_prime, then it should have ulam distance of 0
# BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
# # should be 1
#
#
# test_2 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 2)
# print(test_2)
# # 1 2 3 4 5 6
#
# BayesMallows:::get_rank_distance(rho, test_2, metric= "ulam")
# # should be 0
#
#
# test_3 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 3)
# print(test_3)
# # 1 3 4 2 5 6
#
# BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
# # should be 1
#
#
# # we have a ranking data set containing 10 rankings over 6 items
# test_4 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rankings, metric = metric, rho = rho, leap_size = 1)
# print(test_4)
# # 1 2 3 5 4 6
#
# BayesMallows:::get_rank_distance(rho, test_4, metric= "ulam")
# # should be 1
#
#
#
