#' Estimate Partition Function with Importance Sampling
#'
#' \code{estimate_partition_function} uses importance sampling to estimate
#' the partition function of the Mallows rank model.
#'
#' @param alpha_vector Numeric vector of \eqn{\alpha} values over which
#' to compute the importance sampling estimate.
#'
#' @param n_items Integer specifying the number of items.
#' @param metric Character string specifying the distance measure to use.
#' Available options are \code{"footrule"}, \code{"spearman"},
#' \code{"kendall"}, \code{"cayley"}, and \code{"hamming"}.
#'
#' @param nmc Integer specifying the number of Monte Carlo samples to use in the
#' importance sampling.
#'
#' @param degree Integer specifying the degree of the polynomial used to estimate
#' \eqn{\log(\alpha)} from the grid of values provided by the importance sampling
#' estimate.
#'
#' @return A vector of length \code{degree}.
#'
#' @export
#'
estimate_partition_function <- function(alpha_vector, n_items, metric,
                                        nmc, degree){

  message("This function is not yet implemented.")
}
