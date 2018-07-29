#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesMallows, .registration = TRUE
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
validate_permutation <- function(vec){
  all(sort(vec) == seq_along(vec))
}

#' Create rankings from orderings.
#'
#' @param O A vector or matrix of ordered items. If a matrix,
#' it should be of size N times n, where N is the number of samples
#' and n is the number of items.
#'
#' @return A vector or matrix of rankings.
#' @export
#'
#'
#' @examples
#' # A vector of ordered items.
#' O <- c(5, 1, 2, 4, 3)
#' # Get ranks
#' R <- create_ranking_from_ordering(O)
#' # R is c(2, 3, 5, 4, 1)
#'
#' # Next, we have a matrix with N = 19 samples
#' # and n = 4 items
#' set.seed(21)
#' N <- 10
#' n <- 4
#' O <- t(replicate(N, sample.int(n)))
#' # Convert the ordering to ranking
#' R <- create_ranking_from_ordering(O)
#'
create_ranking_from_ordering <- function(O){

  # Check that it is a permutation
  if(is.vector(O)){
    stopifnot(validate_permutation(O))
    return(order(O))
  } else if(is.matrix(O)){
    # add code for checking matrix
    stopifnot(all(apply(O, 1, validate_permutation)))
    return(t(apply(O, 1, order)))
  } else {
    stop("O must be a vector or matrix")
  }
}
