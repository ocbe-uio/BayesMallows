#' Frequency distribution of the ranking sequences
#'
#' @description Construct the frequency distribution of the distinct ranking sequences from the dataset of the individual rankings.
#'
#' @param rankings A matrix with the individual rankings in each row.
#' @return Numeric matrix with the distinct rankings in each row and the corresponding frequencies indicated in the last \code{(n_items+1)}-th column.
#' @export
#'
#' @example /inst/examples/rank_freq_distr_example.R
#'

rank_freq_distr <- function(rankings){

  if(!is.matrix(rankings)){
    rankings <- matrix(rankings, nrow = 1)
  }

  rankings[which(is.na(rankings))] <- 0
  out <- PLMIX::unit_to_freq(data=rankings)
  out[which(out==0)]=NA

  return(out)

}
