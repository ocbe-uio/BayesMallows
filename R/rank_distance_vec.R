#' Distance between a set of rankings and a given rank sequence
#'
#' @description Compute the distance between a matrix of rankings and a rank sequence.
#' @param rankings A matrix of size \eqn{N }\eqn{\times}{x}\eqn{ n_items} of rankings in each row.
#' @param rho A ranking sequence.
#' @param metric Character string specifying the distance measure to use. Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"}, \code{"ulam"}, \code{"footrule"} and \code{"spearman"}.
#' @return A vector of distances according to the given \code{metric}.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/rank_distance_vec_example.R

rank_distance_vec <- function(rankings, rho, metric){

  rankings <- t(rankings)
  out <- c(rank_dist_vec(rankings=rankings, rho=rho, metric=metric))
  return(out)

}
