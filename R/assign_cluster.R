#' Assign Assessors to Clusters
#'
#' Assign assessors to clusters by finding the cluster with highest
#' posterior probability.
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#'
#' @param soft A logical specifying whether to perform soft or
#' hard clustering. If \code{soft=TRUE}, all cluster probabilities
#' are returned, whereas if \code{soft=FALSE}, only the maximum a
#' posterior (MAP) cluster probability is returned, per assessor. In the
#' case of a tie between two or more cluster assignments, a random cluster
#' is taken as MAP estimate.
#'
#' @param expand A logical specifying whether or not to expand the rowset
#' of each assessor to also include clusters for which the assessor has
#' 0 a posterior assignment probability. Only used when \code{soft = TRUE}. Defaults
#' to \code{FALSE}.
#'
#' @return A dataframe. If \code{soft = FALSE}, it has one row per assessor, and columns \code{assessor},
#' \code{probability} and \code{map_cluster}. If \code{soft = TRUE}, it has \code{n_cluster}
#' rows per assessor, and the additional column \code{cluster}.
#'
#' @seealso \code{\link{compute_mallows}} for an example where this function is used.
#'
#' @export
#'
assign_cluster <- function(model_fit, burnin = model_fit$burnin, soft = TRUE, expand = FALSE){

  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  df <- dplyr::filter(model_fit$cluster_assignment, .data$iteration > burnin)

  # Compute the probability of each iteration
  df <- dplyr::group_by(df, .data$assessor)
  df <- dplyr::mutate(df, probability = 1 / dplyr::n())
  df <- dplyr::ungroup(df)

  # Compute the probability of each cluster, per assessor
  df <- dplyr::group_by(df, .data$assessor, .data$value)
  df <- dplyr::summarise(df, probability = sum(.data$probability))
  df <- dplyr::ungroup(df)
  df <- dplyr::rename(df, cluster = .data$value)

  if(expand){
    df <- tidyr::complete(
      dplyr::group_by(df, .data$assessor),
      cluster = unique(df$cluster),
      fill = list(probability = 0)
    )
    df <- dplyr::ungroup(df)
  }

  # Compute the MAP estimate per assessor
  map <- dplyr::group_by(df, .data$assessor)
  map <- dplyr::mutate(map, max_prob = max(.data$probability))
  map <- dplyr::filter(map, .data$probability == .data$max_prob)

  # Deal with the possible case of ties
  map <- dplyr::filter(map, dplyr::row_number() == 1)
  map <- dplyr::ungroup(map)
  map <- dplyr::select(map, -.data$probability, -.data$max_prob)
  map <- dplyr::rename(map, map_cluster = .data$cluster)

  # Join map back onto df
  df <- dplyr::inner_join(df, map, by = "assessor")

  if(!soft){
    df <- dplyr::filter(df, .data$cluster == .data$map_cluster)
    df <- dplyr::select(df, -.data$cluster)
  }

  return(df)
}
