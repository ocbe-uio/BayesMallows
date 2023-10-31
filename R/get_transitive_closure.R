#' Generic function for obtaining transitive closure
#'
#' @export
#'
get_transitive_closure <- function(rank_data) {
  UseMethod("get_transitive_closure")
}


#' Get transitive closure
#'
#' @param rank_data An object of class \code{"BayesMallowsData"} returned from
#'   \code{\link{setup_rank_data}}.
#'
#' @return A dataframe with transitive closure, if these is any.
#' @export
#'
get_transitive_closure.BayesMallowsData <- function(rank_data) {
  if(inherits(rank_data$preferences, "BayesMallowsTC")) {
    rank_data$preferences
  } else {
    NULL
  }
}
