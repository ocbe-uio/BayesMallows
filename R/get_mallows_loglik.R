#' Likelihood and log-likelihood evaluation for a Mallows mixture model
#'
#' @description Compute either the likelihood or the log-likelihood value of the
#'   Mallows mixture model parameters for a dataset of complete rankings.
#' @param rho A matrix of size \code{n_clusters x n_items} whose rows are
#'   permutations of the first n_items integers corresponding to the modal
#'   rankings of the Mallows mixture components.
#' @param alpha A vector of \code{n_clusters} non-negative scalar specifying the
#'   scale (precision) parameters of the Mallows mixture components.
#' @param weights A vector of \code{n_clusters} non-negative scalars specifying
#'   the mixture weights.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"},
#'   \code{"ulam"} for \code{n_items<=95}, \code{"footrule"} for
#'   \code{n_items<=50} and \code{"spearman"} for \code{n_items<=14}.
#' @param rankings A matrix with observed rankings in each row.
#' @param obs_freq A vector of observation frequencies (weights) to apply to
#'   each row in \code{rankings}. This can speed up computation if a large
#'   number of assessors share the same rank pattern. Defaults to \code{NULL},
#'   which means that each row of \code{rankings} is multiplied by 1. If
#'   provided, \code{obs_freq} must have the same number of elements as there
#'   are rows in \code{rankings}, and \code{rankings} cannot be \code{NULL}.
#' @param log A logical; if TRUE, the log-likelihood value is returned,
#'   otherwise its exponential. Default is \code{TRUE}.
#'
#' @return The likelihood or the log-likelihood value corresponding to one or
#'   more observed complete rankings under the Mallows mixture rank model with
#'   distance specified by the \code{metric} argument.
#' @export
#'
#' @example inst/examples/get_mallows_loglik_example.R
#' @family rank functions
#'
get_mallows_loglik <- function(rho, alpha, weights, metric,
                               rankings, obs_freq = NULL, log = TRUE) {
  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1)
  }

  if (!is.null(obs_freq)) {
    if (nrow(rankings) != length(obs_freq)) {
      stop("obs_freq must be of same length as the number of rows in rankings")
    }
  } else {
    obs_freq <- rep(1, nrow(rankings))
  }

  if (!is.matrix(rho)) {
    rho <- matrix(rho, nrow = 1)
  }

  n_clusters <- length(weights)
  n_items <- ncol(rankings)
  N <- sum(obs_freq)

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
      -(alpha[g] / n_items * rank_dist_sum(
        rankings = t(rankings),
        rho = rho[g, ],
        metric = metric, obs_freq = obs_freq
      ) +
        N * get_partition_function(
          alpha = alpha[g],
          n_items = n_items, metric = metric,
          cardinalities = card
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

#' @rdname get_mallows_loglik
#' @export
lik_db_mix <- function(rho, alpha, weights, metric,
                       rankings, obs_freq = NULL, log = FALSE) {
  .Deprecated(new = "get_mallows_loglik")
  get_mallows_loglik(rho, alpha, weights, metric, rankings, obs_freq, log)
}
