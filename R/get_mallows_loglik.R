#' Likelihood and log-likelihood evaluation for a Mallows mixture model
#'
#' @description Compute either the likelihood or the log-likelihood value of the
#'   Mallows mixture model parameters for a dataset of complete rankings.
#' @param rho A matrix of size `n_clusters x n_items` whose rows are
#'   permutations of the first n_items integers corresponding to the modal
#'   rankings of the Mallows mixture components.
#' @param alpha A vector of `n_clusters` non-negative scalar specifying the
#'   scale (precision) parameters of the Mallows mixture components.
#' @param weights A vector of `n_clusters` non-negative scalars specifying
#'   the mixture weights.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are `"kendall"`, `"cayley"`, `"hamming"`,
#'   `"ulam"`, `"footrule"`, and `"spearman"`.
#' @param rankings A matrix with observed rankings in each row.
#' @param observation_frequency A vector of observation frequencies (weights) to apply to
#'   each row in `rankings`. This can speed up computation if a large
#'   number of assessors share the same rank pattern. Defaults to `NULL`,
#'   which means that each row of `rankings` is multiplied by 1. If
#'   provided, `observation_frequency` must have the same number of elements as there
#'   are rows in `rankings`, and `rankings` cannot be `NULL`.
#' @param log A logical; if TRUE, the log-likelihood value is returned,
#'   otherwise its exponential. Default is `TRUE`.
#'
#' @return The likelihood or the log-likelihood value corresponding to one or
#'   more observed complete rankings under the Mallows mixture rank model with
#'   distance specified by the `metric` argument.
#' @export
#'
#' @example inst/examples/get_mallows_loglik_example.R
#' @family rank functions
#'
get_mallows_loglik <- function(
    rho, alpha, weights,
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
    rankings, observation_frequency = NULL, log = TRUE) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley", "hamming",
    "kendall", "ulam"
  ))
  if (!is.matrix(rankings)) rankings <- matrix(rankings, nrow = 1)
  if (!is.null(observation_frequency)) {
    if (nrow(rankings) != length(observation_frequency)) {
      stop(
        "observation_frequency must be ",
        "of same length as the number of rows in rankings"
      )
    }
  } else {
    observation_frequency <- rep(1, nrow(rankings))
  }

  if (!is.matrix(rho)) rho <- matrix(rho, nrow = 1)

  n_clusters <- length(weights)
  n_items <- ncol(rankings)
  N <- sum(observation_frequency)

  pfun_values <- prepare_partition_function(metric, n_items)

  pp <- sapply(1:n_clusters, function(g) {
    weights[g] * exp(-alpha[g] / n_items * compute_rank_distance(rankings, rho[g, ],
      metric = metric,
      observation_frequency = observation_frequency
    ) - get_partition_function(alpha = alpha[g], n_items = n_items, metric = metric, pfun_values))
  })


  loglik <- sum(log(apply(pp, 1, sum)))


  if (!log) {
    exp(loglik)
  } else {
    loglik
  }
}
