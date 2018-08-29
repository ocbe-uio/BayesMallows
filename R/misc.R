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



gather_rho <- function(model_fit, selected_items = NULL,
                       row_inds = NULL){

  all_items <- seq(1, model_fit$n_items)

  if(!is.null(selected_items)){
    stopifnot(all(selected_items %in% all_items) ||
                all(selected_items %in% rownames(model_fit$rho)))
  } else {
    selected_items <- all_items
  }

  if(!is.null(row_inds)){
    stopifnot(all(row_inds %in% 1:ncol(model_fit$rho)))
  } else {
    row_inds <- seq.int(from = 1, to = ncol(model_fit$rho), by = 1)
  }


  df <- dplyr::tibble(
    index = numeric(),
    cluster = character(),
    item = character(),
    rank = numeric()
  )
  for(i in seq(from = 1, to = model_fit$n_clusters, by = 1)){
    tmp <- matrix(nrow = model_fit$nmc, ncol = length(selected_items))
    tmp[] <- t(model_fit$rho[selected_items,,i])

    if(is.character(selected_items)){
      colnames(tmp) <- selected_items
    } else {
      colnames(tmp) <- rownames(model_fit$rho)[selected_items]
    }


    tmp <- dplyr::as_tibble(tmp)
    tmp <- dplyr::mutate(tmp,
                         index = dplyr::row_number(),
                         cluster = paste("Cluster", i))
    tmp <- tidyr::gather(tmp, key = "item", value = "rank", -.data$index, -.data$cluster)

    df <- dplyr::bind_rows(df, tmp)
  }

  # Convert the item to factor, so it ends up in the order specified
  if(is.character(selected_items)){
    fct_levels <- selected_items
  } else {
    fct_levels <- rownames(model_fit$rho)[selected_items]
  }

  df <- dplyr::mutate(df, item = factor(.data$item, levels = fct_levels))

  return(df)
}


# Function for getting an x axis without decimals.
# Modified from https://stackoverflow.com/questions/21061653/creating-a-density-histogram-in-ggplot2
scalefun <- function(x) sprintf("%d", as.integer(x))


# Function which goes through each row of the tibble pair_comp and
# checks that ranks are valid
validate_initial_ranking <- function(pair_comp, mat){
  pair_comp <- dplyr::rowwise(pair_comp)
  pair_comp <- dplyr::mutate(pair_comp,
                             rank_bottom = mat[[.data$assessor, .data$bottom_item]],
                             rank_top = mat[[.data$assessor, .data$top_item]]
  )

  error_rows <- dplyr::filter(pair_comp, .data$rank_bottom < .data$rank_top)

  if(nrow(error_rows) > 0){
    stop(paste("Invalid ranking generated"))
  } else {
    return(TRUE)
  }
}


prepare_alpha_df <- function(alpha_matrix){
  df <- dplyr::as_tibble(alpha_matrix)
  names(df) <- paste("Cluster", seq(from = 1, to = ncol(alpha_matrix), by = 1))
  df <- dplyr::mutate(df, index = dplyr::row_number())
  df <- tidyr::gather(df, key = "cluster", value = "alpha", -.data$index)

  return(df)
}



#' Validate that rankings are consistent with ordering
#'
#' We need a function to check that ranks are consistent with linear ordering We
#' assume the linear_ordering has n_assessors elements, each being a vector.
#'
#'
#' @param rankings A rank matrix of size n_assessors x n_items.
#'
#' @param linear_ordering A list of length n_assessors. Each list element must
#'   be a vector of orderings. It is assumed that the first vector element is
#'   least favored, and the last element is most preferred.
#'
#' @return A logical. `TRUE` if consistent, `FALSE` otherwise.
#' @export
#'
validate_rank_ordering <- function(rankings, linear_ordering){
  # Convert rankings to list
  rankings <- split(rankings, seq(1, nrow(rankings)))

  # Remove the non-constrained elements
  rankings_constrained <- purrr::map2(linear_ordering, rankings, ~ .y[.x])

  # Verify that all list elements are in ascending order
  all(purrr::map_lgl(rankings_constrained, ~ all(order(.x, decreasing = TRUE) == seq(1, length(.x)))))
}

