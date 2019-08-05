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
#'
#' @return A matrix of rankings which can be given in the \code{rankings} argument
#' to \code{\link{compute_mallows}}.
#'
#' @export
#'
#' @example /inst/examples/generate_initial_ranking_example.R
#'
generate_initial_ranking <- function(tc,
                                     n_items = max(tc[, c("bottom_item", "top_item")]),
                                     cl = NULL){


  if(!("BayesMallowsTC" %in% class(tc))){
    stop("tc must be an object returned from generate_transitive_closure")
  }
  stopifnot(is.null(cl) || inherits(cl, "cluster"))

  prefs <- split(tc[, c("bottom_item", "top_item"), drop = FALSE], tc$assessor)
  if(is.null(cl)){
    prefs <- lapply(prefs, function(x, y) create_ranks(as.matrix(x), y), n_items)
  } else {
    prefs <- parallel::parLapply(cl = cl, X = prefs,
                                 fun = function(x, y) create_ranks(as.matrix(x), y), n_items)
  }

  do.call(rbind, prefs)
}

create_ranks <- function(mat, n_items){
  g <- igraph::graph_from_edgelist(mat)
  g <- as.integer(igraph::topological.sort(g))

  # Add unranked elements at the end
  all_items <- seq(from = 1, to = n_items, by = 1)
  g <- c(g, setdiff(all_items, g))

  # Convert from ordering to ranking
  r <- create_ranking(rev(g))
  mat <- matrix(r, nrow = 1)

  return(mat)
}
