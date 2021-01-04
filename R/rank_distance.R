#' Distance between two rankings
#'
#' @description Compute the distance between two rankings according to several metrics.
#' @details Note that the Spearman distance is the squared L2 norm, whereas
#' the footrule distance is the L1 norm.
#'
#' The Ulam distance uses the SUBSET library developed by John Burkardt, available at http:#'people.sc.fsu.edu/~jburkardt/cpp_src/subset/subset.html.
#'
#' The implementation of Cayley distance is based on a \code{C++} translation of \code{Rankcluster::distCayley} \insertCite{Grimonprez2016}{BayesMallows}.
#' @param r1 A vector of ranks.
#' @param r2 A vector of ranks.
#' @param metric Character string specifying the distance measure to use. Available options are \code{"kendall"}, \code{"cayley"}, \code{"hamming"}, \code{"ulam"}, \code{"footrule"} and \code{"spearman"}.
#' @return A scalar equal to the distance between the two ranking sequences according to the given \code{metric}.
#' @export
#'
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/rank_distance_example.R

rank_distance <- function(r1, r2, metric){

  out <- get_rank_distance(r1=r1, r2=r2, metric=metric)
  return(out)

}
