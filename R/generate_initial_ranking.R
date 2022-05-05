#' Generate Initial Ranking
#'
#' Given a consistent set of pairwise preferences, generate a complete ranking
#' of items which is consistent with the preferences.
#'
#' @param tc A dataframe with pairwise comparisons of \code{S3} subclass
#'   \code{BayesMallowsTC}, returned from
#'   \code{\link{generate_transitive_closure}}.
#'
#' @param n_items The total number of items. If not provided, it is assumed to
#'   equal the largest item index found in \code{tc}, i.e., \code{max(tc[,
#'   c("bottom_item", "top_item")])}.
#'
#' @param cl Optional computing cluster used for parallelization, returned from
#'   \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
#' @param shuffle_unranked Logical specifying whether or not to randomly
#'   permuted unranked items in the initial ranking. When
#'   \code{shuffle_unranked=TRUE} and \code{random=FALSE}, all unranked items
#'   for each assessor are randomly permuted. Otherwise, the first ordering
#'   returned by \code{igraph::topo_sort()} is returned.
#'
#' @param random Logical specifying whether or not to use a random initial
#'   ranking. Defaults to \code{FALSE}. Setting this to \code{TRUE} means that
#'   all possible orderings consistent with the stated pairwise preferences are
#'   generated for each assessor, and one of them is picked at random.
#'
#' @param random_limit Integer specifying the maximum number of items allowed
#'   when all possible orderings are computed, i.e., when \code{random=TRUE}.
#'   Defaults to \code{8L}.
#'
#'
#' @return A matrix of rankings which can be given in the \code{rankings}
#'   argument to \code{\link{compute_mallows}}.
#'
#' @note Setting \code{random=TRUE} means that all possible orderings of each
#'   assessor's preferences are generated, and one of them is picked at random.
#'   This can be useful when experiencing convergence issues, e.g., if the MCMC
#'   algorithm does not mixed properly. However, finding all possible orderings
#'   is a combinatorial problem, which may be computationally very hard. The
#'   result may not even be possible to fit in memory, which may cause the R
#'   session to crash. When using this option, please try to increase the size
#'   of the problem incrementally, by starting with smaller subsets of the
#'   complete data. An example is given below.
#'
#'   As detailed in the documentation to \code{\link{generate_transitive_closure}},
#'   it is assumed that the items are labeled starting from 1. For example, if a single
#'   comparison of the following form is provided, it is assumed that there is a total
#'   of 30 items (\code{n_items=30}), and the initial ranking is a permutation of these 30
#'   items consistent with the preference 29<30.
#'
#' \tabular{rrr}{
#' \strong{assessor} \tab \strong{bottom_item} \tab \strong{top_item}\cr
#' 1 \tab 29 \tab 30\cr
#' }
#'
#' If in reality there are only two items, they should be relabeled to 1 and 2, as follows:
#'
#' \tabular{rrr}{
#' \strong{assessor} \tab \strong{bottom_item} \tab \strong{top_item}\cr
#' 1 \tab 1 \tab 2\cr
#' }
#'
#'
#' @export
#'
#' @example /inst/examples/generate_initial_ranking_example.R
#'
generate_initial_ranking <- function(tc,
                                     n_items = max(tc[, c("bottom_item", "top_item")]),
                                     cl = NULL, shuffle_unranked = FALSE, random = FALSE,
                                     random_limit = 8L) {
  if (!(inherits(tc, "BayesMallowsTC"))) {
    stop("tc must be an object returned from generate_transitive_closure")
  }
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  if (n_items > random_limit && random) {
    stop(paste(
      "Number of items exceeds the limit for generation of random permutations,\n",
      "modify the random_limit argument to override this.\n"
    ))
  }

  if (n_items < max(tc[, c("bottom_item", "top_item")])) {
    stop("Too few items specified. Please see documentation Note about labeling of items.\n")
  }

  prefs <- split(tc[, c("bottom_item", "top_item"), drop = FALSE], tc$assessor)
  if (is.null(cl)) {
    prefs <- lapply(
      prefs, function(x, y, sr, r) create_ranks(as.matrix(x), y, sr, r),
      n_items, shuffle_unranked, random
    )
  } else {
    prefs <- parallel::parLapply(
      cl = cl, X = prefs,
      fun = function(x, y, sr, r) create_ranks(as.matrix(x), y, sr, r),
      n_items, shuffle_unranked, random
    )
  }

  do.call(rbind, prefs)
}

create_ranks <- function(mat, n_items, shuffle_unranked, random) {
  if (!random) {
    g <- igraph::graph_from_edgelist(mat)
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
    r <- create_ranking(rev(g_final))
    mat <- matrix(r, nrow = 1)
  } else {
    graph <- list()
    for (i in seq_len(n_items)) {
      graph[[i]] <- unique(mat[mat[, "top_item"] == i, "bottom_item"])
    }
    indegree_init <- rep(0, n_items)
    indegree <- table(unlist(graph))
    indegree_init[as.integer(names(indegree))] <- indegree
    attr(graph, "indegree") <- indegree_init
    rm(indegree, indegree_init)
    discovered <- rep(FALSE, n_items)
    path <- numeric()

    res <- utils::capture.output(all_topological_sorts(graph, path, discovered, n_items))
    res <- gsub("\\[1\\] ", "", res)
    res <- res[sample(length(res), 1)]

    mat <- matrix(as.numeric(strsplit(res, split = "[^0-9]+")[[1]]), nrow = 1)
  }


  return(mat)
}
