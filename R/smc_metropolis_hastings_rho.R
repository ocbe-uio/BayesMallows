#' @title Metropolis-Hastings Rho
#' @description Function to perform Metropolis-Hastings for new rho under the Mallows model with footrule distance metric!
#' @inheritParams get_mallows_loglik
#' @param leap_size Integer specifying the step size of the leap-and-shift
#' proposal distribution.
#' @export
#' @author Anja Stein
#' @examples
#' rho <- t(c(1,2,3,4,5,6))
#' alpha <- 2
#' metric <- "footrule"
#' n_items <- 6
#'
#' metropolis_hastings_rho(
#' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
#' 	rho = rho, leap_size = 1
#' )
#'
#' metropolis_hastings_rho(
#' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
#' 	rho = rho, leap_size = 2
#' )
#'
#' metropolis_hastings_rho(
#' 	alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
#' 	rho = rho, leap_size = 3
#' )
#'
#' rankings <- sample_mallows(
#'  rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
#' )
#' metropolis_hastings_rho(
#' 	alpha = alpha, n_items = n_items, rankings = rankings, metric = metric,
#' 	rho = rho, leap_size = 1
#' )
#'

metropolis_hastings_rho <- function(alpha, n_items, rankings, metric, rho, leap_size) {


  # create new potential consensus ranking
  kernel <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = leap_size)

  # output from leap-and-shift is of the following
  # leap_shift_list <- list("rho_prime" = rho_prime, "forwards_prob" = forwards_prob, "backwards_prob" = backwards_prob)
  rho_prime <- kernel$rho_prime
  forwards_prob <- kernel$forwards_prob # rho_prime|rho
  backwards_prob <- kernel$backwards_prob # rho|rho_prime

  # evaluate the log-likelihood with current rankings
  mallows_loglik_curr <- get_mallows_loglik(alpha = alpha, rho = rho, n_items = n_items, rankings = rankings, metric = metric)
  mallows_loglik_prop <- get_mallows_loglik(alpha = alpha, rho = rho_prime, n_items = n_items, rankings = rankings, metric = metric)

  # calculate acceptance probability
  loga <- log(backwards_prob) - log(forwards_prob) + mallows_loglik_prop - mallows_loglik_curr


  # determine whether to accept or reject proposed rho and return now consensus ranking
  p <- runif(1, min = 0, max = 1)
  if (log(p) <= loga) {
    return(rho_prime)
  } else {
    return(rho)
  }
}
