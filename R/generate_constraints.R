#' Generate Constraint Set from Pairwise Comparisons
#'
#' This function is relevant when \code{\link{compute_mallows}} is called
#' repeatedly with the same data, e.g., when determining the
#' number of clusters. It precomputes a list of constraints used
#' internally by the MCMC algorithm, which otherwise would be
#' recomputed each time \code{\link{compute_mallows}} is called.
#'
#' @param preferences Data frame of preferences. For the case of consistent
#' rankings, \code{preferences} should be returned from \code{\link{generate_transitive_closure}}.
#' For the case of inconsistent preferences, when using an error model as described
#' in \insertCite{crispino2019;textual}{BayesMallows}, a dataframe of preferences
#' can be directly provided.
#' @param n_items Integer specifying the number of items.
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
#' @return A list which is used internally by the MCMC algorithm.
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/generate_constraints_example.R
#'
generate_constraints <- function(preferences, n_items, cl = NULL){

  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  # Turn the preferences dataframe into a list of dataframes,
  # one list element per assessor
  constraints <- split(preferences[, c("bottom_item", "top_item"), drop = FALSE], preferences$assessor)
  if(is.null(cl)) {
    lapply(constraints, constraint_fun, n_items)
  } else {
    parallel::parLapply(cl = cl, X = constraints, fun = constraint_fun, n_items)
  }
}



constraint_fun <- function(x, n_items){
  # Find out which items are constrained
  constrained_items <- unique(c(x[["bottom_item"]], x[["top_item"]]))

  # Now we must complete the dataframe with the items that do not appear
  items_above <- tidyr::complete(dplyr::select(x, .data$bottom_item, .data$top_item),
                                 bottom_item = seq(from = 1, to = n_items, by = 1),
                                 fill = list(top_item = NA_integer_))

  # Split it into a list, with one element per bottom_item
  items_above <- split(items_above, items_above[["bottom_item"]])

  # For each item, find which items are ranked above it
  items_above <- purrr::map(items_above, function(x) {
    res <- unique(x[["top_item"]])
    res <- res[!is.na(res)]
  })

  # Now we must complete the dataframe with the items that do not appear
  items_below <- tidyr::complete(dplyr::select(x, .data$bottom_item, .data$top_item),
                                 top_item = seq(from = 1, to = n_items, by = 1),
                                 fill = list(bottom_item = NA_integer_))

  # Split it into a list, with one element per bottom_item
  items_below <- split(items_below, items_below[["top_item"]])

  # For each item, find which items are ranked above it
  items_below <- purrr::map(items_below, function(x) {
    res <- unique(x[["bottom_item"]])
    res <- res[!is.na(res)]
  })

  return(
    list(
      constrained_items = constrained_items,
      items_above = items_above,
      items_below = items_below
    )
  )
}
