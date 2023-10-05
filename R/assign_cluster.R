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
#' @family posterior quantities
#'
assign_cluster <- function(model_fit, burnin = model_fit$burnin, soft = TRUE, expand = FALSE) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  if (is.null(model_fit$cluster_assignment)) {
    stop("Rerun compute_mallows with save_clus=TRUE.")
  }
  stopifnot(burnin < model_fit$nmc)

  df <- model_fit$cluster_assignment[model_fit$cluster_assignment$iteration > burnin, , drop = FALSE]

  # Compute the probability of each iteration
  df <- aggregate(
    list(count = df$iteration),
    list(assessor = df$assessor, cluster = df$value),
    FUN = length, drop = !expand
  )
  df$count[is.na(df$count)] <- 0

  df <- do.call(rbind, lapply(split(df, f = df$assessor), function(x) {
    x$probability <- x$count / sum(x$count)
    x$count <- NULL
    x
  }))

  # Compute the MAP estimate per assessor
  map <- do.call(rbind, lapply(split(df, f = df$assessor), function(x) {
    x <- x[x$probability == max(x$probability), , drop = FALSE]
    x <- x[1, , drop = FALSE] # in case of ties
    x$map_cluster <- x$cluster
    x$cluster <- x$probability <- NULL
    x
  }))

  # Join map back onto df
  df <- merge(df, map, by = "assessor")

  if (!soft) {
    df <- df[df$cluster == df$map_cluster, , drop = FALSE]
    df$cluster <- NULL
  }

  return(df)
}
