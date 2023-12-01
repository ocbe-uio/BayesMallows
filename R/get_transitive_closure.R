#' Get transitive closure
#'
#' A simple method for showing any transitive closure computed by
#' [setup_rank_data()].
#'
#' @param rank_data An object of class `"BayesMallowsData"` returned from
#'   [setup_rank_data].
#'
#' @return A dataframe with transitive closure, if there is any.
#' @export
#'
#' @family preprocessing
#'
#' @examples
#' dat <- setup_rank_data(preferences = beach_preferences)
#' head(get_transitive_closure(dat))
#'
get_transitive_closure <- function(rank_data) {
  if (inherits(rank_data$preferences, "BayesMallowsTransitiveClosure")) {
    rank_data$preferences
  } else {
    NULL
  }
}
