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


get_cardinalities <- function(relevant_params, n) {
  unlist(
    dplyr::pull(
      dplyr::select(
        dplyr::filter(
          relevant_params, .data$num_items == n
        ),
        .data$values)
    )
  )
}
