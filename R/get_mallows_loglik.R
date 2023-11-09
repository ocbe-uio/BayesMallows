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
#'   `"ulam"` for `n_items<=95`, `"footrule"` for
#'   `n_items<=50` and `"spearman"` for `n_items<=14`.
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
get_mallows_loglik <- function(rho, alpha, weights, metric,
                               rankings, observation_frequency = NULL, log = TRUE) {
  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1)
  }

  if (!is.null(observation_frequency)) {
    if (nrow(rankings) != length(observation_frequency)) {
      stop("observation_frequency must be of same length as the number of rows in rankings")
    }
  } else {
    observation_frequency <- rep(1, nrow(rankings))
  }

  if (!is.matrix(rho)) {
    rho <- matrix(rho, nrow = 1)
  }

  n_clusters <- length(weights)
  n_items <- ncol(rankings)
  N <- sum(observation_frequency)

  if (metric %in% c("ulam", "footrule", "spearman")) {
    pfd <- partition_function_data[
      partition_function_data$metric == metric &
        partition_function_data$n_items == n_items &
        partition_function_data$type == "cardinalities", ,
      drop = FALSE
    ]
    if (nrow(pfd) == 0) {
      stop("Given number of items currently not available for the specified metric")
    } else {
      card <- pfd$values[[1]]
    }
  } else if (metric %in% c("kendall", "cayley", "hamming")) {
    card <- NULL
  }

  loglik <- vapply(
    X = seq_len(n_clusters),
    FUN = function(g) {
      -(alpha[g] / n_items * sum(rank_dist_vec(
        rankings = t(rankings),
        rho = rho[g, ],
        metric = metric, observation_frequency = observation_frequency
      )) +
        N * get_partition_function(
          alpha = alpha[g],
          n_items = n_items, metric = metric,
          logz_list = list(logz_estimate = NULL, cardinalities = card)
        )) * weights[[g]]
    },
    FUN.VALUE = numeric(1)
  )

  if (!log) {
    exp(sum(loglik))
  } else {
    sum(loglik)
  }
}
