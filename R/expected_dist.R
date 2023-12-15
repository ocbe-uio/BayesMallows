#' Expected value of metrics under a Mallows rank model
#'
#' @description Compute the expectation of several metrics under the Mallows
#'   rank model.
#' @param alpha Non-negative scalar specifying the scale (precision) parameter
#'   in the Mallows rank model.
#' @param n_items Integer specifying the number of items.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are `"kendall"`, `"cayley"`, `"hamming"`, `"ulam"`,
#'   `"footrule"`, and `"spearman"`.
#'
#' @return A scalar providing the expected value of the `metric` under the
#'   Mallows rank model with distance specified by the `metric` argument.
#' @export
#'
#' @family rank functions
#'
#' @example /inst/examples/expected_dist_example.R
compute_expected_distance <- function(
    alpha,
    n_items,
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley",
    "hamming", "kendall", "ulam"
  ))
  validate_integer(n_items)
  validate_positive(n_items)

  if (alpha < 0) {
    stop("alpha must be a non-negative value")
  } else {
    pfun_values <- prepare_partition_function(metric, n_items)
    get_expected_distance(alpha, n_items, metric, pfun_values)
  }
}
