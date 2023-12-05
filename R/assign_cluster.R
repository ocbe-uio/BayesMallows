#' @title Assign Assessors to Clusters
#'
#' @description Assign assessors to clusters by finding the cluster with highest
#' posterior probability.
#'
#' @param model_fit An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param burnin A numeric value specifying the number of iterations to discard
#'   as burn-in. Defaults to `model_fit$burnin`, and must be provided if
#'   `model_fit$burnin` does not exist. See [assess_convergence()].
#'
#' @param soft A logical specifying whether to perform soft or hard clustering.
#'   If `soft=TRUE`, all cluster probabilities are returned, whereas if
#'   `soft=FALSE`, only the maximum a posterior (MAP) cluster probability is
#'   returned, per assessor. In the case of a tie between two or more cluster
#'   assignments, a random cluster is taken as MAP estimate.
#'
#' @param expand A logical specifying whether or not to expand the rowset of
#'   each assessor to also include clusters for which the assessor has 0 a
#'   posterior assignment probability. Only used when `soft = TRUE`. Defaults to
#'   `FALSE`.
#'
#' @return A dataframe. If `soft = FALSE`, it has one row per assessor, and
#'   columns `assessor`, `probability` and `map_cluster`. If `soft = TRUE`, it
#'   has `n_cluster` rows per assessor, and the additional column `cluster`.
#'
#' @export
#'
#' @family posterior quantities
#'
#' @examples
#' # Fit a model with three clusters to the simulated example data
#' set.seed(1)
#' mixture_model <- compute_mallows(
#'   data = setup_rank_data(cluster_data),
#'   model_options = set_model_options(n_clusters = 3),
#'   compute_options = set_compute_options(nmc = 5000, burnin = 1000)
#' )
#'
#' head(assign_cluster(mixture_model))
#' head(assign_cluster(mixture_model, soft = FALSE))
#'
assign_cluster <- function(
    model_fit, burnin = model_fit$burnin, soft = TRUE, expand = FALSE) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  df <- model_fit$cluster_assignment[
    model_fit$cluster_assignment$iteration > burnin, , drop = FALSE]

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

  map <- do.call(rbind, lapply(split(df, f = df$assessor), function(x) {
    x <- x[x$probability == max(x$probability), , drop = FALSE]
    x <- x[1, , drop = FALSE] # in case of ties
    x$map_cluster <- x$cluster
    x$cluster <- x$probability <- NULL
    x
  }))

  df <- merge(df, map, by = "assessor")

  if (!soft) {
    df <- df[df$cluster == df$map_cluster, , drop = FALSE]
    df$cluster <- NULL
  }

  return(df)
}
