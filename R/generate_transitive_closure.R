#' Generate Transitive Closure
#'
#' Generate the transitive closure for a set of consistent pairwise comparisons. The result
#' can be given in the \code{preferences} argument to \code{\link{compute_mallows}}.
#'
#' @param df A data frame with one row per pairwise comparison, and columns
#' \code{assessor}, \code{top_item}, and \code{bottom_item}. Each column contains the
#' following:
#' \itemize{
#' \item \code{assessor} is a numeric vector containing the assessor index, or a character
#'       vector containing the (unique) name of the assessor.
#'
#' \item \code{bottom_item} is a numeric vector containing the index of the item that
#'       was disfavored in each pairwise comparison.
#'
#' \item \code{top_item} is a numeric vector containing the index of the item that was
#'       preferred in each pairwise comparison.
#' }
#' So if we have two assessors and five items, and assessor 1 prefers item 1 to item 2 and
#' item 1 to item 5, while assessor 2 prefers item 3 to item 5, we have the following \code{df}:
#' \tabular{rrr}{
#' \strong{assessor} \tab \strong{bottom_item} \tab \strong{top_item}\cr
#' 1 \tab 2 \tab 1\cr
#' 1 \tab 5 \tab 1\cr
#' 2 \tab 5 \tab 3\cr
#' }
#'
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}.
#'
#'
#' @return A dataframe with the same columns as \code{df}, but with its set of rows expanded
#' to include all pairwise preferences implied by the ones stated in \code{df}. The returned
#' object has \code{S3} subclass \code{BayesMallowsTC}, to indicate that this is the
#' transitive closure.
#'
#' @seealso \code{\link{generate_initial_ranking}}
#'
#' @export
#'
#' @example /inst/examples/generate_transitive_closure_example.R
#'
generate_transitive_closure <- function(df, cl = NULL){

  stopifnot(is.null(cl) || inherits(cl, "cluster"))
  prefs <- split(df[, c("bottom_item", "top_item"), drop = FALSE], df$assessor)

  if(is.null(cl)){
    prefs <- mapply(function(x, y) cbind(y, .generate_transitive_closure(cbind(x$bottom_item, x$top_item))),
                    prefs, unique(df$assessor), SIMPLIFY = FALSE)
  } else {
    prefs <- parallel::clusterMap(cl = cl,
                                  fun = function(x, y) cbind(y, .generate_transitive_closure(cbind(x$bottom_item, x$top_item))),
                                  prefs, unique(df$assessor))
  }


  prefs <- do.call(rbind.data.frame, prefs)
  rownames(prefs) <- NULL
  colnames(prefs) <- colnames(df)

  # Check if there are any inconsistencies
  check <- dplyr::semi_join(prefs, prefs,
                            by = c("assessor" = "assessor",
                                   "bottom_item" = "top_item",
                                   "top_item" = "bottom_item"))

  if(nrow(check) > 0){
    print("Inconsistent rankings:")
    print(check)
    stop("Cannot compute transitive closure. Please run compute_mallows with error_model='bernoulli'.")
  }

  class(prefs) <- c("BayesMallowsTC", class(prefs))

  return(prefs)
}


#' Internal function for generating transitive closure
#'
#' @param mat A matrix in which column 1 is the lower ranked item and column 2 is the
#'   upper ranked item.
#' @keywords internal
.generate_transitive_closure <- function(mat){

  # This line was an answer to StackOverflow question 51794127
  my_set <- do.call(sets::set, apply(mat, 1, sets::as.tuple))

  # Next I compute the transitive closure:
  r <- relations::endorelation(graph = my_set)
  tc <- relations::transitive_closure(r)
  incidence <- relations::relation_incidence(tc)

  new_mat <- which(incidence == 1, arr.ind = TRUE)

  result <- data.frame(
    bottom_item = as.numeric(rownames(incidence)[new_mat[, 1, drop = FALSE]]),
    top_item = as.numeric(colnames(incidence)[new_mat[, 2, drop = FALSE]])
    )


  return(result)

}

