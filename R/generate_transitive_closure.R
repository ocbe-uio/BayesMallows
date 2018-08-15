#' Generate transitive closure
#'
#' @param df Data frame with columns \code{assessor}, \code{top_item}, and \code{bottom_item}.
#' @export
generate_transitive_closure <- function(df){

  # Integers mess up the function, so we use numeric
  df <- dplyr::mutate_all(df, as.numeric)

  df <- dplyr::group_by(df, .data$assessor)
  result <- dplyr::do(
    df,
    dplyr::as_tibble(
      .generate_transitive_closure(
        cbind(.data$bottom_item, .data$top_item))
      )
  )
  result <- dplyr::ungroup(result)

  names(result) <- names(df)
  class(result) <- c("BayesMallowsTC", class(result))

  return(result)
}


#' Internal function for generating transitive closure
#'
#' @param mat A matrix in which column 1 is the lower ranked item and column 2 is the
#'   upper ranked item.
.generate_transitive_closure <- function(mat){

  # This line was an answer to StackOverflow question 51794127
  my_set <- do.call(sets::set, apply(mat, 1, sets::as.tuple))

  # Next I compute the transitive closure:
  r <- relations::endorelation(graph = my_set)
  tc <- relations::transitive_closure(r)
  incidence <- relations::relation_incidence(tc)

  new_mat <- which(incidence == 1, arr.ind = TRUE)

  result <- cbind(
    as.integer(rownames(incidence)[new_mat[, 1, drop = FALSE]]),
    as.integer(colnames(incidence)[new_mat[, 2, drop = FALSE]])
  )

  return(result)

}
