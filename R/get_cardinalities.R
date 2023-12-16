#' @title Get cardinalities for each distance
#'
#' @description The partition function for the Mallows model can be defined in a
#'   computationally efficient manner as
#'   \deqn{Z_{n}(\alpha) = \sum_{d_{n} \in
#'   \mathcal{D}_{n}} N_{m,n} e^{-(\alpha/n) d_{m}}}.
#'   In this equation, \eqn{\mathcal{D}_{n}} a set containing all possible
#'   distances at the given number of items, and \eqn{d_{m}} is on element of
#'   this set. Finally, \eqn{N_{m,n}} is the number of possible configurations
#'   of the items that give the particular distance. See
#'   \insertCite{irurozki2016;textual}{BayesMallows},
#'   \insertCite{vitelli2018;textual}{BayesMallows}, and
#'   \insertCite{crispino2023;textual}{BayesMallows} for details.
#'
#'   For footrule distance, the cardinalities come from entry A062869 in the
#'   On-Line Encyclopedia of Integer Sequences (OEIS)
#'   \insertCite{oeis}{BayesMallows}. For Spearman distance, they come from
#'   entry A175929, and for Ulam distance from entry A126065.
#'
#' @param n_items Number of items.
#' @param metric Distance function, one of "footrule", "spearman", or "ulam".
#'
#' @return A dataframe with two columns, `distance` which contains each distance
#'   in the support set at the current number of items, i.e., \eqn{d_{m}}, and
#'   `value` which contains the number of values at this particular distances,
#'   i.e., \eqn{N_{m,n}}.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example inst/examples/get_cardinalities_example.R
#' @family partition function
get_cardinalities <- function(
    n_items,
    metric = c("footrule", "spearman", "ulam")) {
  metric <- match.arg(metric, c("footrule", "spearman", "ulam"))

  if (metric == "footrule") {
    if (n_items > length(footrule_cardinalities)) {
      stop("Not available for requested number of items.")
    }
    as.data.frame(footrule_cardinalities[[n_items]])
  } else if (metric == "spearman") {
    if (n_items > length(spearman_cardinalities)) {
      stop("Not available for requested number of items.")
    }
    as.data.frame(spearman_cardinalities[[n_items]])
  } else if (metric == "ulam") {
    if (n_items > length(ulam_cardinalities)) {
      stop("Not available for requested number of items.")
    }
    as.data.frame(ulam_cardinalities[[n_items]])
  }
}

#' @title Distances and cardinalities
#'
#' @description List containing distances and cardinalities for footrule
#'   distance. See [get_cardinalities()] for details.
#'
#' @family partition function
"footrule_cardinalities"

#' @title Distances and cardinalities
#'
#' @description List containing distances and cardinalities for Spearman
#'   distance. See [get_cardinalities()] for details.
#'
#' @family partition function
"spearman_cardinalities"

#' @title Distances and cardinalities
#'
#' @description List containing distances and cardinalities for Ulam
#'   distance. See [get_cardinalities()] for details.
#'
#' @family partition function
"ulam_cardinalities"
