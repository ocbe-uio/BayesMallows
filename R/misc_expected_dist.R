exp_d_tau <- function(alpha, n_items) {
  if (alpha > 0) {
    idx <- seq(from = 1, to = n_items, by = 1)
    out <- n_items * exp(-alpha) / (1 - exp(-alpha)) -
      sum((idx * exp(-idx * alpha)) / (1 - exp(-idx * alpha)))
  } else {
    if (alpha == 0) {
      out <- n_items * (n_items - 1) / 4
    } else {
      stop("alpha must be a non-negative value")
    }
  }
  return(out)
}

exp_d_cay <- function(alpha, n_items) {
  idx <- seq(from = 1, to = n_items - 1, by = 1)
  out <- sum(idx / (idx + exp(alpha)))
  return(out)
}

exp_d_ham <- function(alpha, n_items) {
  idx <- seq(from = 0, to = n_items, by = 1)
  out <- n_items - exp(alpha) *
    sum(((exp(alpha) - 1)^idx[-(n_items + 1)]) / factorial(idx[-(n_items + 1)])) /
    sum(((exp(alpha) - 1)^idx) / factorial(idx))
  return(out)
}

