splitpref <- function(preferences) {
  lapply(split(
    preferences[, c("bottom_item", "top_item"), drop = FALSE],
    preferences$assessor
  ), as.matrix)
}

generate_initial_ranking <- function(
    preferences, n_items, cl = NULL, max_topological_sorts) {
  UseMethod("generate_initial_ranking")
}

#' @export
generate_initial_ranking.BayesMallowsTransitiveClosure <- function(
    preferences, n_items, cl = NULL, max_topological_sorts) {
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  prefs <- splitpref(preferences)

  if (is.null(cl)) {
    do.call(rbind, lapply(
      prefs, function(x, y, r) create_ranks(x, y, r),
      n_items, max_topological_sorts
    ))
  } else {
    do.call(rbind, parallel::parLapply(
      cl = cl, X = prefs,
      fun = function(x, y, r) create_ranks(x, y, r),
      n_items, max_topological_sorts
    ))
  }
}

#' @export
generate_initial_ranking.BayesMallowsIntransitive <- function(
    preferences, n_items, cl = NULL, max_topological_sorts) {
  n_assessors <- length(unique(preferences$assessor))
  rankings <- replicate(n_assessors, sample(x = n_items, size = n_items),
    simplify = "numeric"
  )
  rankings <- matrix(rankings, ncol = n_items, nrow = n_assessors, byrow = TRUE)
}

create_ranks <- function(mat, n_items, max_topological_sorts) {
  ret <- all_topological_sorts(mat, n_items, max_topological_sorts, TRUE)
  u <- sample(min(max_topological_sorts, nrow(ret)), 1)
  ret <- ret[u, ]
  all_items <- seq(from = 1, to = n_items, by = 1)
  ranked_items <- unique(c(mat))
  unranked_items <- setdiff(all_items, ranked_items)
  idx_ranked <- sort(sample(length(all_items), length(ranked_items)))
  g_final <- rep(NA, n_items)
  g_final[idx_ranked] <- ret[ret %in% ranked_items]
  g_final[is.na(g_final)] <- unranked_items[sample(length(unranked_items))]
  create_ranking(g_final)
}
