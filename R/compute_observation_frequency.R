#' Frequency distribution of the ranking sequences
#'
#' @description Construct the frequency distribution of the distinct ranking
#'   sequences from the dataset of the individual rankings. This can be of
#'   interest in itself, but also used to speed up computation by providing
#'   the `observation_frequency` argument to [compute_mallows()].
#'
#' @param rankings A matrix with the individual rankings in each row.
#' @return Numeric matrix with the distinct rankings in each row and the
#'   corresponding frequencies indicated in the last `(n_items+1)`-th
#'   column.
#' @export
#' @family rank functions
#'
#' @example /inst/examples/compute_observation_frequency_example.R
#'
compute_observation_frequency <- function(rankings) {
  if(!is.matrix(rankings)) stop("rankings must be a matrix")

  rankings[is.na(rankings)] <- 0
  counts <- table(apply(rankings, 1, paste, collapse = ","))

  ret <- cbind(
    do.call(rbind,
            lapply(strsplit(names(counts), split = ","), as.numeric)),
    as.numeric(counts)
  )
  ret[ret == 0] <- NA
  ret

}

