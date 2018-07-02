#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesMallows, .registration = TRUE
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @examples
#' n <- 10
#' vec <- sample(n,n)
#' validate_permutation(vec)
validate_permutation <- function(vec){
  all(sort(vec) == seq_along(vec))
}
