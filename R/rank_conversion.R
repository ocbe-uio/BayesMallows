#' Convert between ranking and ordering.
#'
#' \code{create_ranking} takes a vector or matrix of ordered items \code{orderings} and
#' returns a corresponding vector or matrix of ranked items.
#' \code{create_ordering} takes a vector or matrix of rankings \code{rankings} and
#' returns a corresponding vector or matrix of ordered items.
#'
#' @param orderings A vector or matrix of ordered items. If a matrix, it should be of
#'   size N times n, where N is the number of samples and n is the number of
#'   items.
#' @param rankings A vector or matrix of ranked items. If a matrix, it should be N
#'   times n, where N is the number of samples and n is the number of items.
#'
#' @return A vector or matrix of rankings. Missing orderings coded as \code{NA} are propagated into corresponding missing ranks and vice versa.
#'
#'
#' @examples
#' # A vector of ordered items.
#' orderings <- c(5, 1, 2, 4, 3)
#' # Get ranks
#' rankings <- create_ranking(orderings)
#' # rankings is c(2, 3, 5, 4, 1)
#' # Finally we convert it backed to an ordering.
#' orderings_2 <- create_ordering(rankings)
#' # Confirm that we get back what we had
#' all.equal(orderings, orderings_2)
#'
#' # Next, we have a matrix with N = 19 samples
#' # and n = 4 items
#' set.seed(21)
#' N <- 10
#' n <- 4
#' orderings <- t(replicate(N, sample.int(n)))
#' # Convert the ordering to ranking
#' rankings <- create_ranking(orderings)
#' # Now we try to convert it back to an ordering.
#' orderings_2 <- create_ordering(rankings)
#' # Confirm that we get back what we had
#' all.equal(orderings, orderings_2)
#' @name rank_conversion
NULL

#' @describeIn rank_conversion Convert from ordering to ranking.
#' @export
create_ranking <- function(orderings){

  # Check that it is a permutation
  if(is.vector(orderings)){
    stopifnot(validate_permutation(orderings))
    return(order(orderings))
  } else if(is.matrix(orderings)){
    n_items <- ncol(orderings)

    # Convert to list, for easier functional programming
    orderings <- purrr::array_branch(orderings, margin = 1)

    # Check that matrix contains permutations
    check <- purrr::map_lgl(orderings, validate_permutation)
    if(!all(check)){
      stop(paste("orderings must contain proper permutations. Problem row(s):", which(!check)))
    }

    # Convert each ordering to ranking, taking special care of missing values
    rankings <- purrr::map(orderings, function(x) {
      # Find out which items are missing
      missing_items <- setdiff(1:n_items, x)
      # Possible rankings for each item
      candidates <- matrix(1:n_items, ncol = n_items, nrow = n_items)
      candidates[, missing_items] <- NA_real_
      # Logical matrix specifying which item to pick
      inds <- outer(X = x, Y = 1:n_items, FUN = "==")
      inds[1, colSums(inds, na.rm = TRUE) == 0] <- TRUE
      inds[is.na(inds)] <- FALSE
      # Extract the correct items
      candidates[inds]
    } )

    return(t(matrix(unlist(rankings), ncol = length(rankings))))
  } else {
    stop("orderings must be a vector or matrix")
  }
}

#' @describeIn rank_conversion Convert from ranking to ordering.
#' @export
create_ordering <- function(rankings){

  create_ranking(rankings)
}
