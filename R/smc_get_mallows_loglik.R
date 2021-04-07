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
#' @param rankings # TODO: describe
#' @return Mallows log-likelihood
#' @export
#' @examples
#' set.seed(101)
#' rho <- c(1,2,3,4,5,6)
#' alpha <- 2
#' metric <- "footrule"
#' n_items <- 6
#' get_mallows_loglik(
#'   alpha = alpha, rho = rho, n_items = length(rho), rankings = rho,
#'   metric = metric
#' )
#'
#' # return 0 because you are comparing the consensus ranking with itself
#' # if you change alpha or metric, then the result shall remain as 0
#'
#' rankings <- sample_mallows(
#'   rho0 = rho, alpha0 = alpha, n_samples = 10, burnin = 1000, thinning = 500
#' )
#'
#' # depending on your seed, you will get a different collection of rankings in R and C++
#'
#' get_mallows_loglik(
#'   alpha = alpha, rho = rho,  n_items = n_items, rankings = rankings ,
#'   metric = metric
#' )
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
