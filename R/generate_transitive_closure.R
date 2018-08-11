#' Title
#'
#' @param df Dataframe with columns \code{Assessor}, \code{BottomItem}, and
#'   \code{TopItem}.
#' @param num_items Number of items in total.
#'
#' @export
#'
generate_transitive_closure <- function(df, num_items){

  df <- dplyr::group_by(df, .data$Assessor)
  result <- dplyr::do(
    df,
    dplyr::as_tibble(
      .generate_transitive_closure(
        cbind(.data$BottomItem, .data$TopItem),
        num_items = num_items)
      )
  )
  result <- dplyr::ungroup(result)

  names(result) <- names(df)

  # Add the number of items as an attribute
  attr(result, "num_items") <- num_items

  return(result)
}


#' Internal function for generating transitive closure
#'
#' @param mat A matrix in which column 1 is the lower ranked item and column 2 is the
#'   upper ranked item.
#' @param num_items Total number of items.
#'
#'
#' @return
#'
#' @examples
.generate_transitive_closure <- function(mat, num_items){

  incidence_matrix <- matrix(0, nrow = num_items, ncol = num_items)
  colnames(incidence_matrix) <- seq(from = 1, to = num_items, by = 1)
  rownames(incidence_matrix) <- seq(from = 1, to = num_items, by = 1)

  # This line was an answer to StackOverflow question 51794127
  my_set <- do.call(sets::set, apply(mat, 1, sets::as.tuple))

  # Next I compute the transitive closure:
  r <- relations::endorelation(graph = my_set)
  tc <- relations::transitive_closure(r)
  incidence <- relations::relation_incidence(tc)

  # Put this into the incidence_matrix defined above
  incidence_matrix[rownames(incidence), colnames(incidence)] <- incidence

  # We now have a new matrix, with potentially more rows that the original
  result <- arrayInd(which(incidence_matrix == 1), .dim = c(num_items, num_items))

  return(result)

}
