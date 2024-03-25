#' @title Compute exact partition function
#'
#' @description For Cayley, Hamming, and Kendall distances, computationally
#' tractable functions are available for the exact partition function.
#'
#' @param alpha Dispersion parameter.
#' @param n_items Number of items.
#' @param metric Distance function, one of "cayley", "hamming", or "kendall".
#'
#' @return The logarithm of the partition function.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example inst/examples/compute_exact_partition_function_example.R
#' @family partition function
compute_exact_partition_function <- function(
    alpha, n_items,
    metric = c("cayley", "hamming", "kendall")) {
  metric <- match.arg(metric, c("cayley", "hamming", "kendall"))
  validate_integer(n_items)
  validate_positive(n_items)
  validate_positive(alpha)

  get_partition_function(alpha, n_items, metric, NULL)
}
