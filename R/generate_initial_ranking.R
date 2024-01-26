generate_initial_ranking <- function(
    preferences, cl = NULL, shuffle_unranked = FALSE,
    random = FALSE, random_limit = 8L) {
  UseMethod("generate_initial_ranking")
}

#' @export
generate_initial_ranking.BayesMallowsTransitiveClosure <- function(
    preferences, cl = NULL, shuffle_unranked = FALSE, random = FALSE,
    random_limit = 8L) {
  stopifnot(is.null(cl) || inherits(cl, "cluster"))
  n_items <- max(preferences[, c("bottom_item", "top_item")])
  if (n_items > random_limit && random) {
    stop(paste(
      "Number of items exceeds the limit for generation of random permutations,\n",
      "modify the random_limit argument to override this.\n"
    ))
  }

  prefs <- split(
    preferences[, c("bottom_item", "top_item"), drop = FALSE],
    preferences$assessor
  )
  if (is.null(cl)) {
    do.call(rbind, lapply(
      prefs, function(x, y, sr, r) create_ranks(as.matrix(x), y, sr, r),
      n_items, shuffle_unranked, random
    ))
  } else {
    do.call(rbind, parallel::parLapply(
      cl = cl, X = prefs,
      fun = function(x, y, sr, r) create_ranks(as.matrix(x), y, sr, r),
      n_items, shuffle_unranked, random
    ))
  }
}

#' @export
generate_initial_ranking.BayesMallowsIntransitive <- function(
    preferences, cl = NULL, shuffle_unranked = FALSE,
    random = FALSE, random_limit = 8L) {
  n_items <- max(preferences[, c("bottom_item", "top_item")])
  n_assessors <- length(unique(preferences$assessor))
  rankings <- replicate(n_assessors, sample(x = n_items, size = n_items),
    simplify = "numeric"
  )
  rankings <- matrix(rankings, ncol = n_items, nrow = n_assessors, byrow = TRUE)
}

.S3method(
  "generate_initial_ranking", "BayesMallowsIntransitive",
  generate_initial_ranking.BayesMallowsIntransitive
)

create_ranks <- function(mat, n_items, shuffle_unranked, random) {
  if (!random) {
    g <- igraph::graph_from_edgelist(as.matrix(mat))
    g <- as.integer(igraph::topo_sort(g))

    all_items <- seq(from = 1, to = n_items, by = 1)

    if (!shuffle_unranked) {
      # Add unranked elements outside of the range at the end
      g_final <- c(g, setdiff(all_items, g))
    } else {
      ranked_items <- unique(c(mat))
      unranked_items <- setdiff(all_items, ranked_items)
      # Indices of ranked elements in final vector
      idx_ranked <- sort(sample(length(all_items), length(ranked_items)))
      g_final <- rep(NA, n_items)
      g_final[idx_ranked] <- g[g %in% ranked_items]
      g_final[is.na(g_final)] <- unranked_items[sample(length(unranked_items))]
    }

    # Convert from ordering to ranking
    return(create_ranking(rev(g_final)))
  } else {
    graph <- list()
    for (i in seq_len(n_items)) {
      graph[[i]] <- unique(mat[mat[, "top_item"] == i, "bottom_item"])
    }
    indegree_init <- rep(0, n_items)
    indegree <- table(unlist(graph))
    indegree_init[as.integer(names(indegree))] <- indegree
    attr(graph, "indegree") <- indegree_init

    e1 <- new.env()
    assign("x", list(), envir = e1)
    assign("num", 0L, envir = e1)
    all_topological_sorts(graph, n_items, e1)
    return(get("x", envir = e1)[[sample(get("num", envir = e1), 1)]])
  }
}
