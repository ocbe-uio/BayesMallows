.onUnload <- function(libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @noRd
validate_permutation <- function(vec) {
  if (!any(is.na(vec))) {
    return(all(sort(vec) == seq_along(vec)))
  } else if (all(is.na(vec))) {
    return(TRUE)
  } else {
    return(all(vec[!is.na(vec)] <= length(vec)) &&
      all(vec[!is.na(vec)] >= 1) && !any(duplicated(vec[!is.na(vec)])))
  }
}

## Source: https://stackoverflow.com/questions/11095992/generating-all-distinct-permutations-of-a-list-in-r
permutations <- function(n) {
  if (n == 1) {
    return(matrix(1))
  } else {
    sp <- permutations(n - 1)
    p <- nrow(sp)
    A <- matrix(nrow = n * p, ncol = n)
    for (i in 1:n) {
      A[(i - 1) * p + 1:p, ] <- cbind(i, sp + (sp >= i))
    }
    return(A)
  }
}
