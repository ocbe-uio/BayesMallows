splitpref <- function(preferences) {
  lapply(split(
    preferences[, c("bottom_item", "top_item"), drop = FALSE],
    preferences$assessor
  ), as.matrix)
}

generate_initial_ranking <- function(
    preferences, n_items, cl = NULL, random = FALSE, random_limit = 8L) {
  UseMethod("generate_initial_ranking")
}

#' @export
generate_initial_ranking.BayesMallowsTransitiveClosure <- function(
    preferences, n_items, cl = NULL, random = FALSE,
    random_limit = 8L) {
  stopifnot(is.null(cl) || inherits(cl, "cluster"))
  if (n_items > random_limit && random) {
    stop(paste(
      "Number of items exceeds the limit for generation of random permutations,\n",
      "modify the random_limit argument to override this.\n"
    ))
  }

  prefs <- splitpref(preferences)

  if (is.null(cl)) {
    do.call(rbind, lapply(
      prefs, function(x, y, r) create_ranks(x, y, r),
      n_items, random
    ))
  } else {
    do.call(rbind, parallel::parLapply(
      cl = cl, X = prefs,
      fun = function(x, y, r) create_ranks(x, y, r),
      n_items, random
    ))
  }
}

#' @export
generate_initial_ranking.BayesMallowsIntransitive <- function(
    preferences, n_items, cl = NULL,
    random = FALSE, random_limit = 8L) {
  n_assessors <- length(unique(preferences$assessor))
  rankings <- replicate(n_assessors, sample(x = n_items, size = n_items),
    simplify = "numeric"
  )
  rankings <- matrix(rankings, ncol = n_items, nrow = n_assessors, byrow = TRUE)
}

create_ranks <- function(mat, n_items, random) {
  if (!random) {
    g <- igraph::graph_from_edgelist(mat)
    g <- as.integer(igraph::topo_sort(g))

    all_items <- seq(from = 1, to = n_items, by = 1)
    ranked_items <- unique(c(mat))
    unranked_items <- setdiff(all_items, ranked_items)
    # Indices of ranked elements in final vector
    idx_ranked <- sort(sample(length(all_items), length(ranked_items)))
    g_final <- rep(NA, n_items)
    g_final[idx_ranked] <- g[g %in% ranked_items]
    g_final[is.na(g_final)] <- unranked_items[sample(length(unranked_items))]

    # Convert from ordering to ranking
    return(create_ranking(rev(g_final)))
  } else {
    ret <- all_topological_sorts(mat, n_items)
    create_ranking(ret[sample(nrow(ret), 1), ])
  }
}
