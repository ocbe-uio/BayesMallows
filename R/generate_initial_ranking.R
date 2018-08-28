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
#' @return A matrix of rankings which can be given in the \code{rankings} argument
#' to \code{\link{compute_mallows}}.
#'
#' @export
#'
#' @example /inst/examples/generate_initial_ranking_example.R
#'
generate_initial_ranking <- function(tc,
                                     n_items = max(tc[, c("bottom_item", "top_item")])){

  if(!("BayesMallowsTC" %in% class(tc))){
    stop("tc must be an object returned from generate_transitive_closure")
  }

  # Find
  tc <- dplyr::group_by(tc, .data$assessor)
  tc <- dplyr::do(tc, dplyr::as_tibble(
      x = create_ranks(
        as.matrix(.data[, c("bottom_item", "top_item"), drop = FALSE]),
        n_items = n_items)
      ))

  mat <- as.matrix(tc[, -1, drop = FALSE])
  rownames(mat) <- tc[["assessor"]]
  colnames(mat) <- seq(from = 1, to = n_items, by = 1)
  return(mat)
}

create_ranks <- function(mat, n_items){
  g <- create_linear_ordering(mat)

  # Add unranked elements at the end
  all_items <- seq(from = 1, to = n_items, by = 1)
  g <- c(g, setdiff(all_items, g))

  # Convert from ordering to ranking
  r <- create_ranking(rev(g))
  mat <- matrix(r, nrow = 1)

  return(mat)
}


create_linear_ordering <- function(mat, partial = FALSE){
  g <- igraph::graph_from_edgelist(mat)
  g <- as.integer(igraph::topological.sort(g))

  if(partial){
    g <- base::intersect(g, c(mat))
  }
  return(g)
}
