#' Distance between a set of rankings and a given rank sequence
#'
#' @description Compute the distance between a matrix of rankings and a rank
#'   sequence.
#' @param rankings A matrix of size \eqn{N \times n_{items}} of
#'   rankings in each row. Alternatively, if \eqn{N} equals 1, `rankings`
#'   can be a vector.
#' @param rho A ranking sequence.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are `"kendall"`, `"cayley"`, `"hamming"`,
#'   `"ulam"`, `"footrule"` and `"spearman"`.
#' @param observation_frequency Vector of observation frequencies of length \eqn{N}, or of length 1,
#' which means that all ranks are given the same weight. Defaults to 1.
#' @return A vector of distances according to the given `metric`.
#' @export
#'
#' @details The implementation of Cayley distance is based on a `C++`
#' translation of `Rankcluster::distCayley()` \insertCite{Grimonprez2016}{BayesMallows}.
#'
#' @references \insertAllCited
#' @family rank functions
#'
#' @example /inst/examples/compute_rank_distance_example.R
compute_rank_distance <- function(
    rankings, rho,
    metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
    observation_frequency = 1) {
  metric <- match.arg(metric, c(
    "footrule", "spearman", "cayley", "hamming",
    "kendall", "ulam"
  ))
  if (!is.matrix(rankings)) rankings <- matrix(rankings, nrow = 1)

  stopifnot(length(observation_frequency) == 1 ||
    length(observation_frequency) == nrow(rankings))

  if (length(observation_frequency) == 1) {
    observation_frequency <- rep(observation_frequency, nrow(rankings))
  }

  as.numeric(
    get_rank_distance(rankings = t(rankings), rho = rho, metric = metric) *
      observation_frequency)
}
