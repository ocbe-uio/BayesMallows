#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesMallows, .registration = TRUE
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @export
#' @examples
#' n <- 10
#' vec <- sample(n,n)
#' isPerm(vec)
validate_permutation <- function(vec){
  n <- length(vec)

  if( (sum(vec) == sum(1:n) ) && ( sum(vec^2) == sum((1:n)^2))){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
