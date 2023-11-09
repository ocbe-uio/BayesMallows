generate_constraints <- function(data, cl = NULL) {
  if (is.null(data$preferences)) {
    return(list())
  }
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  # Turn the preferences dataframe into a list of dataframes,
  # one list element per assessor
  constraints <- split(
    data$preferences[, c("bottom_item", "top_item"), drop = FALSE],
    data$preferences$assessor
  )
  if (is.null(cl)) {
    lapply(constraints, constraint_fun, data$n_items)
  } else {
    parallel::parLapply(cl = cl, X = constraints, fun = constraint_fun, data$n_items)
  }
}



constraint_fun <- function(x, n_items) {
  # Find out which items are constrained
  constrained_items <- unique(c(x[["bottom_item"]], x[["top_item"]]))

  # Now we must complete the dataframe with the items that do not appear
  items_above <- merge(x[, c("bottom_item", "top_item"), drop = FALSE],
    expand.grid(bottom_item = seq(from = 1, to = n_items, by = 1)),
    by = "bottom_item", all = TRUE
  )

  # Split it into a list, with one element per bottom_item
  items_above <- split(items_above, items_above$bottom_item)

  # For each item, find which items are ranked above it
  items_above <- lapply(items_above, function(x) {
    res <- unique(x[["top_item"]])
    res <- res[!is.na(res)]
  })

  # Now we must complete the dataframe with the items that do not appear
  items_below <- merge(x[, c("bottom_item", "top_item"), drop = FALSE],
    expand.grid(top_item = seq(from = 1, to = n_items, by = 1)),
    by = "top_item", all = TRUE
  )

  # Split it into a list, with one element per bottom_item
  items_below <- split(items_below, items_below$top_item)

  # For each item, find which items are ranked above it
  items_below <- lapply(items_below, function(x) {
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
