metropolis_hastings_rho_R <- function(alpha, n_items, rankings, metric, rho, leap_size) {


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
