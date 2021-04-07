#' @title Get Mallows log-likelihood
#' @description Calculates the Mallows log-likelihood given a set of rankings and a given rank sequence
#' @param alpha Numeric value of the scale parameter
#' @param rho A ranking sequence
#' @param n_items Integer is the number of items in a ranking
#' A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
#' rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
#' can be a vector.
#' @param metric Character string specifying the distance measure to use.
#' Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"},
#' \code{"ulam"}, \code{"footrule"} and \code{"spearman"}.
#' @return Mallows log-likelihood
#' @export
get_mallows_loglik <- function(alpha, rho, n_items, rankings, metric) {
  sum_distance <- 0
  num_rankings <- dim(rankings)[1]

  # calculate the sum of the distances
  if (is.null(num_rankings)) {
    sum_distance <- sum_distance + get_rank_distance(rho, rankings, metric = metric)
  } else {
    for (jj in 1:num_rankings) {
      sum_distance <- sum_distance + get_rank_distance(rho, rankings[jj, ], metric = metric)
    }
  }

  mallows_loglik <- -alpha / n_items * sum_distance
  return(mallows_loglik)
}
