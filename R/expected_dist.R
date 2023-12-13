#' Expected value of metrics under a Mallows rank model
#'
#' @description Compute the expectation of several metrics under the Mallows
#'   rank model.
#' @param alpha Non-negative scalar specifying the scale (precision) parameter
#'   in the Mallows rank model.
#' @param n_items Integer specifying the number of items.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are `"kendall"`, `"cayley"`, `"hamming"`, `"ulam"` for
#'   `n_items<=95`, `"footrule"`, and `"spearman"`.
#'
#' @return A scalar providing the expected value of the `metric` under the
#'   Mallows rank model with distance specified by the `metric` argument.
#' @export
#'
#' @family rank functions
#'
#' @example /inst/examples/expected_dist_example.R
compute_expected_distance <- function(
    alpha,
    n_items,
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley",
    "hamming", "kendall", "ulam"
  ))
  if (n_items < 1 | floor(n_items) != n_items) {
    stop("Number of items must be a positive integer")
  }

  # Scale alpha to parametrization used
  alpha <- alpha / n_items

  if (alpha < 0) {
    stop("alpha must be a non-negative value")
  } else {
    if (metric == "kendall") {
      out <- exp_d_tau(alpha, n_items)
    }
    if (metric == "cayley") {
      out <- exp_d_cay(alpha, n_items)
    }
    if (metric == "hamming") {
      out <- exp_d_ham(alpha, n_items)
    }
    if (metric %in% c("ulam", "footrule", "spearman")) {
      card <- get_cardinalities(n_items, metric)$value
      out <- exp(
        log_expected_dist(
          alpha = alpha * n_items,
          n_items = n_items,
          cardinalities = card,
          metric = metric
        )
      )
    }
  }
  return(out)
}

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
