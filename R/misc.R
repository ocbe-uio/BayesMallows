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
