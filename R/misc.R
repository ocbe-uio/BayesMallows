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




prepare_partition_function <- function(logz_estimate, metric, n_items){
  # First, has the user supplied an estimate?
  if(!is.null(logz_estimate)){
    return(list(cardinalities = NULL, logz_estimate = logz_estimate))
  }

  # Second, do we have a sequence?
  relevant_params <- dplyr::filter(partition_function_data, .data$n_items == !!n_items,
                                   .data$metric == !!metric, .data$type == "cardinalities")
  if(nrow(relevant_params) == 1){
    return(list(cardinalities = unlist(relevant_params$values), logz_estimate = NULL))
  }

  # Third, do we have an importance sampling estimate?
  relevant_params <- dplyr::filter(partition_function_data, .data$n_items == !!n_items,
                                   .data$metric == !!metric, .data$type == "importance_sampling")

  if(nrow(relevant_params) == 1){
    return(list(cardinalities = NULL, logz_estimate = unlist(relevant_params$values)))
  }

  # Fourth, is it the Ulam distance?
  if(metric == "ulam"){
    return(list(
      cardinalities = unlist(lapply(0:(n_items - 1),
                                    function(x) PerMallows::count.perms(perm.length = n_items, dist.value = x, dist.name = "ulam")))
    ))
  }

  # Fifth, can we compute the partition function in our C++ code?
  if(metric %in% c("cayley", "hamming", "kendall")){
    return(list(cardinalities = NULL, logz_estimate = NULL))
  }

  stop("Partition function not available. Please compute an estimate using estimate_partition_function().")

}
