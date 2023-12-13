#' @title Get cardinalities for each distance
#'
#' @description The partition function for the Mallows model can be defined in a
#'   computationally efficient manner as
#'   \deqn{Z_{n}(\alpha) = \sum_{d_{n} \in
#'   \mathcal{D}_{n}} N_{m,n} e^{-(\alpha/n) d_{m}}}.
#'   In this equation, \eqn{\mathcal{D}_{n}} a set containing all possible
#'   distances at the given number of items, and \eqn{d_{m}} is on element of
#'   this set. Finally, \eqn{N_{m,n}} is the number of possible
#'   configurations of the items that give the particular distance. See
#'   \insertCite{irurozki2016;textual}{BayesMallows},
#'   \insertCite{vitelli2018;textual}{BayesMallows}, and
#'   \insertCite{crispino2023;textual}{BayesMallows} for details.
#'
#'
#' @param n_items Number of items.
#' @param metric Distance function.
#'
#' @return A dataframe with two columns, `distance` which contains each
#'   distance in the support set at the current number of items, i.e.,
#'   \eqn{d_{m}}, and `value` which contains the number of values at this
#'   particular distances, i.e., \eqn{N_{m,n}}.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @examples
#' # Extract the cardinalities for four items with footrule distance
#' n_items <- 4
#' dat <- get_cardinalities(n_items)
#' # Compute the partition function at alpha = 2
#' alpha <- 2
#' sum(dat$value * exp(-alpha / n_items * dat$distance))
#' @family preprocessing
get_cardinalities <- function(n_items, metric = "footrule") {
  metric <- match.arg(metric, "footrule")

  if(metric == "footrule") {
    if(n_items > length(footrule_cardinalities)) {
      stop("Not available for requested number of items.")
    }
    as.data.frame(footrule_cardinalities[[n_items]])
  }
}
