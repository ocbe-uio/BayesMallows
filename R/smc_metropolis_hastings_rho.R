#' @title Metropolis-Hastings Rho
#' @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!
#'   @param alpha Numeric value og the scale parameter
#'   @param n_items Integer is the number of items in a ranking
#'   @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
#'   rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
#'   can be a vector.
#'   @param rho A ranking sequence vector
#'   @param leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution.
#' @return \code{rho} or \code{rho_prime}: A ranking sequence vector to be the next value of rho in the MCMC chain
metropolis_hastings_rho <- function(alpha, n_items, rankings, metric, rho, leap_size){


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
