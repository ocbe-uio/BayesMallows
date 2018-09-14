#' Generate Constraint Set from Pairwise Comparisons
#'
#' This function is relevant when compute_mallows is called
#' repeatedly with the same data, e.g., when determining the
#' number of clusters. It precomputes a list of constraints used
#' internally by the MCMC algorithm, which otherwise would be
#' recomputed each time compute_mallows is called.
#'
#' @param preferences Data frame of preferences returned from
#' \code{\link{generate_transitive_closure}}.
#' @param n_items Integer specifying the number of items.
#'
#' @return A list which is used internally by the MCMC algorithm.
#' @export
#'
#' @examples
#' # Here is an example with the beach preference data.
#' # First, generate the transitive closure.
#' beach_tc <- generate_transitive_closure(beach_preferences)
#' # Next, generate an initial ranking.
#' beach_init_rank <- generate_initial_ranking(beach_tc)
#' # Then generate the constrain set used intervally by compute_mallows
#' constr <- generate_constraints(beach_tc, n_items = 15)
#' # Provide all this elements to compute_mallows
#' model_fit <- compute_mallows(rankings = beach_init_rank,
#' preferences = beach_tc, constraints = constr)
#'
generate_constraints <- function(preferences, n_items){
  # Turn the preferences dataframe into a list of dataframes,
  # one list element per assessor
  constraints <- split(preferences, preferences$assessor)

  purrr::map(constraints, function(x){
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
  })
}
