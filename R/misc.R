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



gather_rho <- function(model_fit, selected_items = NULL,
                       row_inds = NULL){

  all_items <- seq(1, model_fit$n_items)

  if(!is.null(selected_items)){
    stopifnot(all(selected_items %in% all_items) )
  } else {
    selected_items <- all_items
  }

  if(!is.null(row_inds)){
    stopifnot(all(row_inds %in% 1:ncol(model_fit$rho)))
  } else {
    row_inds <- seq.int(from = 1, to = ncol(model_fit$rho), by = 1)
  }

  df <- dplyr::as_tibble(t(model_fit$rho[selected_items, row_inds, drop = FALSE]))
  names(df) <- selected_items

  df <- dplyr::mutate(df, Index = dplyr::row_number())

  # Make the tibble tall by gathering items
  df <- tidyr::gather(df, key = "Item", value = "Rank", -.data$Index)

  # Convert the item to factor
  df <- dplyr::mutate(df,
                      Item = factor(paste("Item", .data$Item),
                                    levels = paste("Item", sort(selected_items))))

  return(df)
}
