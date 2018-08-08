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


get_cardinalities <- function(relevant_params) {
  unlist(
    dplyr::pull(
      dplyr::select(
        dplyr::filter(
          relevant_params
        ),
        .data$values)
    )
  )
}


gather_rho <- function(model_fit, selected_items = NULL){

  all_items <- seq(1, model_fit$n_items)

  if(!is.null(selected_items)){
    stopifnot(all(selected_items %in% all_items) )
  } else {
    selected_items <- all_items
  }

  df <- dplyr::as_tibble(t(model_fit$rho[selected_items, , drop = FALSE]))
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
