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

  # Initialize an empty set
  my_set <- sets::set()

  # Loop over the rows of df, adding pairs to my_set
  for(i in seq(from = 1, to = nrow(mat), by = 1)){
    my_set <- sets::set_union(my_set, sets::pair(mat[[i, 1]], mat[[i, 2]]))
  }

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

#
# # I have sets of pairwise comparisons of the form "Item 2 is better then Item
# # 1". I want to generate the transitive closure of these pairwise comparisons.
#
# # Here is my current approach, for the example case:
# # {Item1 < Item2, Item2 < Item5, Item4 < Item5}
#
# # Total number of items
# #num_items <- 5
#
# rownames(incidence_matrix) <- colnames(incidence_matrix) <- seq(from = 1, to = num_items, by = 1)
#
# # Generate a tibble to hold the data
# df <- dplyr::tribble(
#   ~'Left', ~'Right',
#   1, 2,
#   2, 5,
#   4, 5
# )
#
# # Initialize an empty set
# my_set <- sets::set()
#
# # Loop over the rows of df, adding pairs to my_set
# for(i in seq(from = 1, to = nrow(df), by = 1)){
#   my_set <- sets::set_union(my_set, sets::pair(df[[i, 1]], df[[i, 2]]))
# }
#
# # I have now converted the tibble to a set
# my_set
#
# # Next I compute the transitive closure:
# r <- relations::endorelation(graph = my_set)
# tc <- relations::transitive_closure(r)
# incidence <- relations::relation_incidence(tc)
# incidence
#
# # Put this into the incidence_matrix defined above
# incidence_matrix[rownames(incidence), colnames(incidence)] <- incidence
# incidence_matrix
#
# # We now have a new matrix, with potentially more rows that the original
# # dataframe df
# full_set <- dplyr::as_tibble(
#   arrayInd(
#     which(incidence_matrix == 1),
#     .dim = c(num_items, num_items)))
#
# names(full_set) <- names(df)
#
# # This is the resulting dataframe
# full_set
