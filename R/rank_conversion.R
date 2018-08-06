#' Convert between ranking and ordering.
#'
#' \code{create_ranking} takes a vector or matrix of ordered items \code{O} and
#' returns a corresponding vector or matrix of ranked items.
#' \code{create_ordering} takes a vector or matrix of rankings \code{R} and
#' returns a corresponding vector or matrix of ordered items.
#'
#' @param O A vector or matrix of ordered items. If a matrix, it should be of
#'   size N times n, where N is the number of samples and n is the number of
#'   items.
#' @param R A vector or matrix of ranked items. If a matrix, it should be N
#'   times n, where N is the number of samples and n is the number of items.
#'
#' @return A vector or matrix of rankings.
#'
#'
#' @examples
#' # A vector of ordered items.
#' O <- c(5, 1, 2, 4, 3)
#' # Get ranks
#' R <- create_ranking(O)
#' # R is c(2, 3, 5, 4, 1)
#' # Finally we convert it backed to an ordering.
#' O2 <- create_ordering(R)
#' # Confirm that we get back what we had
#' all.equal(O, O2)
#'
#' # Next, we have a matrix with N = 19 samples
#' # and n = 4 items
#' set.seed(21)
#' N <- 10
#' n <- 4
#' O <- t(replicate(N, sample.int(n)))
#' # Convert the ordering to ranking
#' R <- create_ranking(O)
#' # Now we try to convert it back to an ordering.
#' O2 <- create_ordering(R)
#' # Confirm that we get back what we had
#' all.equal(O, O2)
#' @name rank_conversion
NULL

#' @describeIn rank_conversion Convert from ordering to ranking.
#' @export
create_ranking <- function(O){

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

#' @describeIn rank_conversion Convert from ranking to ordering.
#' @export
create_ordering <- function(R){
  create_ranking(R)
}
