#' Distance between a set of rankings and a given rank sequence
#'
#' @description Compute the distance between a matrix of rankings and a rank
#'   sequence.
#' @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of
#'   rankings in each row. Alternatively, if \eqn{N} equals 1, \code{rankings}
#'   can be a vector.
#' @param rho A ranking sequence.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"},
#'   \code{"ulam"}, \code{"footrule"} and \code{"spearman"}.
#' @param obs_freq Vector of observation frequencies of length \eqn{N}, or of length 1,
#' which means that all ranks are given the same weight. Defaults to 1.
#' @return A vector of distances according to the given \code{metric}.
#' @export
#'
#' @details The implementation of Cayley distance is based on a \code{C++}
#' translation of \code{Rankcluster::distCayley} \insertCite{Grimonprez2016}{BayesMallows}.
#'
#' @references \insertAllCited
#'
#' @example /inst/examples/rank_distance_example.R
rank_distance <- function(rankings, rho, metric, obs_freq = 1) {

  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1)
  }

  stopifnot(length(obs_freq) == 1 || length(obs_freq) == nrow(rankings))

  if (length(obs_freq) == 1) {
    obs_freq <- rep(obs_freq, nrow(rankings))
  }

  rankings <- t(rankings)
  out <- c(rank_dist_vec(rankings = rankings, rho = rho,
                       metric = metric, obs_freq = obs_freq))
  return(out)

}
