#' Generate Transitive Closure
#'
#' Generate the transitive closure for a set of consistent pairwise comparisons. The result
#' can be given in the `preferences` argument to [compute_mallows()].
#'
#' @param df A data frame with one row per pairwise comparison, and columns
#' `assessor`, `top_item`, and `bottom_item`. Each column contains the
#' following:
#' \itemize{
#' \item `assessor` is a numeric vector containing the assessor index, or a character
#'       vector containing the (unique) name of the assessor.
#'
#' \item `bottom_item` is a numeric vector containing the index of the item that
#'       was disfavored in each pairwise comparison.
#'
#' \item `top_item` is a numeric vector containing the index of the item that was
#'       preferred in each pairwise comparison.
#' }
#' So if we have two assessors and five items, and assessor 1 prefers item 1 to item 2 and
#' item 1 to item 5, while assessor 2 prefers item 3 to item 5, we have the following `df`:
#' \tabular{rrr}{
#' **assessor** \tab **bottom_item** \tab **top_item**\cr
#' 1 \tab 2 \tab 1\cr
#' 1 \tab 5 \tab 1\cr
#' 2 \tab 5 \tab 3\cr
#' }
#'
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from [parallel::makeCluster()]. Defaults to `NULL`.
#'
#'
#' @return A dataframe with the same columns as `df`, but with its set of rows expanded
#' to include all pairwise preferences implied by the ones stated in `df`. The returned
#' object has `S3` subclass `BayesMallowsTransitiveClosure`, to indicate that this is the
#' transitive closure.
#'
#' @seealso [generate_initial_ranking()]
#'
#' @export
#'
#' @example /inst/examples/generate_transitive_closure_example.R
#'
#' @family preprocessing
#'
generate_transitive_closure <- function(df, cl = NULL) {
  if(is.null(df)) return(NULL)
  stopifnot(is.null(cl) || inherits(cl, "cluster"))
  prefs <- split(df[, c("bottom_item", "top_item"), drop = FALSE], df$assessor)

  fun <- function(x, y) cbind(y, .generate_transitive_closure(cbind(x$bottom_item, x$top_item)))
  if (is.null(cl)) {
    prefs <- Map(fun, prefs, unique(df$assessor))
  } else {
    prefs <- parallel::clusterMap(cl = cl, fun = fun, prefs, unique(df$assessor))
  }


  prefs <- do.call(rbind.data.frame, prefs)
  rownames(prefs) <- NULL
  colnames(prefs) <- colnames(df)

  # Check if there are any inconsistencies
  check <- merge(prefs, prefs,
    by.x = c("assessor", "bottom_item", "top_item"),
    by.y = c("assessor", "top_item", "bottom_item")
  )

  if (nrow(check) > 0) {
    stop(
      "Cannot compute transitive closure, but that is fine. Just make ",
      "sure you run compute_mallows with argument error_model='bernoulli'")
  }

  class(prefs) <- c("BayesMallowsTransitiveClosure", class(prefs))

  return(prefs)
}


#' Internal function for generating transitive closure
#'
#' @param mat A matrix in which column 1 is the lower ranked item and column 2 is the
#'   upper ranked item.
#' @noRd
.generate_transitive_closure <- function(mat) {
  # This line was an answer to StackOverflow question 51794127
  my_set <- do.call(sets::set, apply(mat, 1, sets::as.tuple))

  # Next I compute the transitive closure:
  r <- relations::endorelation(graph = my_set)
  tc <- relations::transitive_closure(r)
  incidence <- relations::relation_incidence(tc)

  new_mat <- which(incidence == 1, arr.ind = TRUE)
  row_inds <- as.numeric(gsub("[^0-9]+", "", rownames(incidence)))
  result <- data.frame(
    bottom_item = row_inds[new_mat[, 1, drop = FALSE]],
    top_item = row_inds[new_mat[, 2, drop = FALSE]]
  )


  return(result)
}
