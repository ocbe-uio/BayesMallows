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
