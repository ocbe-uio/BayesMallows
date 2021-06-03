#' Generate Initial Ranking
#'
#' Given a consistent set of pairwise preferences, generate a complete ranking
#' of items which is consistent with the preferences.
#'
#' @param tc A dataframe with pairwise comparisons of \code{S3} subclass
#' \code{BayesMallowsTC}, returned from \code{\link{generate_transitive_closure}}.
#'
#' @param n_items The total number of items. If not provided, it is assumed to
#'   equal the largest item index found in \code{tc}, i.e.,
#'   \code{max(tc[, c("bottom_item", "top_item")])}.
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
#' @param random Logical specifying whether or not to use a random initial ranking.
#'   Defaults to \code{FALSE}.
#'
#'
#' @return A matrix of rankings which can be given in the \code{rankings} argument
#' to \code{\link{compute_mallows}}.
#'
#' @note Setting \code{random=TRUE} means that all possible orderings of each assessor's
#' preferences are generated, and one of them is picked at random. This can be useful
#' when experiencing convergence issues, e.g., if the MCMC algorithm does not mixed
#' properly. However, finding all possible orderings is a combinatorial problem,
#' which may be computationally very hard. The result may not even be possible to fit in
#' memory, which may cause the R session to crash. When using this option,
#' please try to increase the size of the problem incrementally, by starting with smaller
#' subsets of the complete data. An example is given below.
#'
#' @export
#'
#' @example /inst/examples/generate_initial_ranking_example.R
#'
generate_initial_ranking <- function(tc,
                                     n_items = max(tc[, c("bottom_item", "top_item")]),
                                     cl = NULL, random = FALSE){


  if(!("BayesMallowsTC" %in% class(tc))){
    stop("tc must be an object returned from generate_transitive_closure")
  }
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  prefs <- split(tc[, c("bottom_item", "top_item"), drop = FALSE], tc$assessor)
  if(is.null(cl)){
    prefs <- lapply(prefs, function(x, y, r) create_ranks(as.matrix(x), y, r), n_items, random)
  } else {
    prefs <- parallel::parLapply(cl = cl, X = prefs,
                                 fun = function(x, y, r) create_ranks(as.matrix(x), y, r), n_items, random)
  }

  do.call(rbind, prefs)
}

create_ranks <- function(mat, n_items, random){

  if(!random){
    g <- igraph::graph_from_edgelist(mat)
    g <- as.integer(igraph::topo_sort(g))

    # Add unranked elements at the end
    all_items <- seq(from = 1, to = n_items, by = 1)
    g <- c(g, setdiff(all_items, g))

    # Convert from ordering to ranking
    r <- create_ranking(rev(g))
    mat <- matrix(r, nrow = 1)
  } else {
    graph <- list()
    for(i in seq_len(n_items)){
      graph[[i]] <- unique(mat[mat[, "top_item"] == i, "bottom_item"])
    }
    indegree_init <- rep(0, n_items)
    indegree <- table(unlist(graph))
    indegree_init[as.integer(names(indegree))] <- indegree
    attr(graph, "indegree") <- indegree_init
    rm(indegree, indegree_init)
    discovered <- rep(FALSE, n_items)
    path <- numeric()

    stdout <- vector('character')
    con <- textConnection('stdout', 'wr', local = TRUE)
    sink(con)
    all_topological_sorts(graph, path, discovered, n_items)
    sink()
    close(con)

    res <- gsub("\\[1\\] ", "", stdout[sample(length(stdout), 1)])

    mat <- matrix(as.numeric(strsplit(res, split = "[^0-9]+")[[1]]), nrow = 1)
  }


  return(mat)
}
