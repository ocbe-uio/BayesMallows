#' Frequency distribution of the ranking sequences
#'
#' @description Construct the frequency distribution of the distinct ranking
#'   sequences from the dataset of the individual rankings. This can be of
#'   interest in itself, but also used to speed up computation by providing
#'   the `obs_freq` argument to [compute_mallows()].
#'
#' @param rankings A matrix with the individual rankings in each row.
#' @return Numeric matrix with the distinct rankings in each row and the
#'   corresponding frequencies indicated in the last `(n_items+1)`-th
#'   column.
#' @export
#' @family rank functions
#'
#' @example /inst/examples/rank_freq_distr_example.R
#'
rank_freq_distr <- function(rankings) {
  if (!is.matrix(rankings)) {
    rankings <- matrix(rankings, nrow = 1)
  }

  rankings[is.na(rankings)] <- 0
  out <- unit_to_freq(data = rankings)
  out[out == 0] <- NA

  return(out)
}
