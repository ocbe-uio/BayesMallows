#' @importFrom Rdpack reprompt
#' @importFrom Rcpp sourceCpp
#' @useDynLib BayesMallows, .registration = TRUE
NULL


.onUnload <- function (libpath) {
  library.dynam.unload("BayesMallows", libpath)
}

#' Check if a vector is a permutation
#' @param vec a vector
#' @return TRUE if vec is a permutation
#' @keywords internal
validate_permutation <- function(vec){
  if(!any(is.na(vec))){
    return(all(sort(vec) == seq_along(vec)))
  } else if(all(is.na(vec))){
    return(TRUE)
  } else {
    return(all(vec[!is.na(vec)] <= length(vec)) &&
             all(vec[!is.na(vec)] >= 1) && !any(duplicated(vec[!is.na(vec)])))
  }
}




# Function for getting an x axis without decimals.
# Modified from https://stackoverflow.com/questions/21061653/creating-a-density-histogram-in-ggplot2
scalefun <- function(x) sprintf("%d", as.integer(x))




