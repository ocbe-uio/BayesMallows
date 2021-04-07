#' @title Leap and Shift Probabilities
#'  @description Determine the new Calculates transition probabilities for proposing a new rho
#'   @param rho A ranking sequence
#'   @param leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution.
#'   @param n_items Integer is the number of items in a ranking

#' @return A list containing:
#' \itemize{
#'   \item \code{rho_prime} A ranking sequence proposed consensus ranking
#'   \item \code{forwards_prob} Numeric value to account for transition probability from rho to rho_prime
#'   \item \code{backwards_prob} Numeric Value to account for the transition probability from \code{rho_prime} to \code{rho}
#' }
leap_and_shift_probs <- function(rho, leap_size, n_items) {


  # draw u uniformly from {1,...,n}
  u <- sample.int(n_items, 1)

  # define set of integers S - the support set for sampling new ranks
  low_bd <- max(1, (rho[u] - leap_size))
  max_bd <- min(n_items, (rho[u] + leap_size))
  S <- seq(from = low_bd, to = max_bd, by = 1)
  S <- S[S != rho[u]]

  # print(S)
  # draw a random number r from S
  r <- sample(S, 1)

  # Create leap step
  rho_star <- rho
  rho_star[u] <- r # replace uth entry with r


  # here, two elements are the same so we need to shift element and replace the repeated r with u
  delta <- rho_star[u] - rho[u]
  rho_prime <- c(rep(0, n_items))

  # shift step
  for (i in 1:n_items) {
    if (rho[i] == rho[u]) {
      rho_prime[i] <- rho_star[u]
    } else if ((rho[u] < rho[i]) & (rho[i] <= rho_star[u]) & (delta > 0)) {
      rho_prime[i] <- rho[i] - 1
    } else if ((rho[u] > rho[i]) & (rho[i] >= rho_star[u]) & (delta < 0)) {
      rho_prime[i] <- rho[i] + 1
    } else {
      rho_prime[i] <- rho[i]
    }
  }

  # Define support set for ranks rho_star[u] can leap to
  S_star <- seq(from = max(1, (rho_star[u] - leap_size)), to = min(n_items, (rho_star[u] + leap_size)))
  S_star <- S_star[S_star != rho_star[u]]

  # calculate forward and backwards probabilities
  if (abs(rho_star[u] - rho[u]) == 1) {

    # p(proposed|current)
    forwards_prob <- 1 / (n_items * length(S)) + 1 / (n_items * length(S_star))

    # p(current|proposed)
    backwards_prob <- forwards_prob
  } else {

    # p(proposed|current)
    forwards_prob <- 1 / (n_items * length(S))

    # p(current|proposed)
    backwards_prob <- 1 / (n_items * length(S_star))
  }


  # check - rho_prime has a ulam distance of 1 from rho
  # BayesMallows:::get_rank_distance(rho, rho_prime, metric="ulam")

  leap_shift_list <- list("rho_prime" = rho_prime, "forwards_prob" = forwards_prob, "backwards_prob" = backwards_prob)

  return(leap_shift_list)
}
