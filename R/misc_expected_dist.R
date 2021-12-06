#/' exp_d_tau
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Kendall distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Kendall metric under the Mallows rank model with the Kendall distance.

exp_d_tau <- function(alpha, n_items) {
  if(alpha > 0) {
    idx <- seq(from = 1, to = n_items, by = 1)
    out <- n_items * exp(-alpha) / (1 - exp(-alpha)) -
      sum((idx * exp(-idx * alpha)) / (1 - exp(-idx * alpha)))
  } else {
    if(alpha == 0) {
      out <- n_items * (n_items - 1) / 4
    } else {
      stop("alpha must be a non-negative value")
    }
  }
  return(out)
}

#/' exp_d_cay
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Cayley distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Cayley metric under the Mallows rank model with the Cayley distance.

exp_d_cay <- function(alpha, n_items) {
  idx <- seq(from = 1, to = n_items - 1, by = 1)
  out <- sum(idx / (idx + exp(alpha)))
  return(out)
}

#/' exp_d_ham
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Hamming distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Hamming metric under the Mallows rank model with the Hamming distance.

exp_d_ham <- function(alpha, n_items) {
  idx <- seq(from = 0, to = n_items, by = 1)
  out <- n_items - exp(alpha) *
    sum(((exp(alpha) - 1)^idx[-(n_items + 1)]) / base::factorial(idx[-(n_items + 1)])) /
    sum(((exp(alpha) - 1)^idx) / base::factorial(idx))
  return(out)
}

#/' exp_d_ulam
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Ulam distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Ulam metric under the Mallows rank model with the Ulam distance.

# The function based on the command from the PerMallows package is slow because it has to generate the distance frequencies.
exp_d_ulam <- function(alpha, n_items) { # for n_items<=95
  idx <- seq(from = 0, to = n_items - 1, by = 1)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric == "ulam"][[n_items]]
  norm_const <- exp(
    get_partition_function(
      alpha = alpha * n_items, n_items = n_items, metric = "ulam",
      cardinalities = card
    )
  )
  out <- sum(idx * exp(-alpha * idx) * card) / norm_const
  return(out)
}

#/' exp_d_foot
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Footrule distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Footrule metric under the Mallows rank model with the Footrule distance.

exp_d_foot <- function(alpha, n_items) { # for n_items<=50
  idx <- seq(0, floor(n_items^2 / 2), by = 2)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric == "footrule"][[n_items]]
  norm_const <- exp(
    get_partition_function(
      alpha = alpha * n_items, n_items = n_items, metric = "footrule",
      cardinalities = card
    )
  )
  out <- sum(idx * exp(-alpha * idx) * card) / norm_const
  return(out)
}

#/' exp_d_spear
#/'
#/' @param alpha Non-negative scalar specifying the scale (precision) parameter in the Mallows rank model with the Spearman distance.
#/' @param n_items Integer specifying the number of items.
#/'
#/' @return Expected value of the Spearman metric under the Mallows rank model with the Spearman distance.

exp_d_spear <- function(alpha, n_items) { # for n_items<=14
  idx <- seq(0, 2 * base::choose(n_items + 1, 3), by = 2)
  pfd <- partition_function_data
  card <- pfd$values[pfd$metric == "spearman"][[n_items]]
  norm_const <- exp(
    get_partition_function(
      alpha = alpha * n_items, n_items = n_items, metric = "spearman",
      cardinalities = card
    )
  )
  out <- sum(idx * exp(-alpha * idx) * card) / norm_const
  return(out)
}
