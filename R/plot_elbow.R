#' Plot Within-Cluster Sum of Distances
#'
#' Plot the within-cluster sum of distances from the corresponding cluster
#' consensus for different number of clusters. This function is useful for
#' selecting the number of mixture.
#'
#' @param ... One or more objects returned from \code{\link{compute_mallows}},
#' separated by comma, or a list of such objects. Typically, each object
#' has been run with a different number of mixtures, as specified in the
#' \code{n_clusters} argument to \code{\link{compute_mallows}}.
#'
#' @param burnin The number of iterations to discard as burnin. Either a vector of
#' numbers, one for each model, or a single number which is taken to be the burnin for
#' all models. If each model provided has a \code{burnin} element, then this is taken
#' as the default.
#'
#' @return A boxplot with the number of clusters on the horizontal axis and the
#' with-cluster sum of distances on the vertical axis.
#'
#' @export
#'
#' @seealso \code{\link{compute_mallows}}
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#'
plot_elbow <- function(..., burnin = NULL) {
  # Put the models into a list. These are typically fitted with different number of clusters
  models <- list(...)

  # Taking into account the case where the user has entered a list of models
  if (length(models) == 1 && !(inherits(models[[1]], "BayesMallows"))) {
    models <- models[[1]]
  }

  df <- do.call(rbind, lapply(models, function(x) {
    stopifnot(inherits(x, "BayesMallows"))

    if (!("burnin" %in% names(x))) {
      if (is.null(burnin)) {
        stop("burnin not provided")
      } else {
        x$burnin <- burnin
      }
    }

    if (!x$include_wcd) stop("To get an elbow plot, set include_wcd=TRUE in compute_mallows")

    df <- x$within_cluster_distance[x$within_cluster_distance$iteration > x$burnin, , drop = FALSE]

    if (nrow(df) <= 0) stop("burnin must be strictly smaller than the number of MCMC samples")

    # Need to sum the within-cluster distances across clusters, for each iteration
    df <- aggregate(x = list(value = df$value), by = list(iteration = df$iteration), FUN = sum)

    df$n_clusters <- x$n_clusters
    return(df)
  }))

  ggplot2::ggplot(df, ggplot2::aes(x = as.factor(.data$n_clusters), y = .data$value)) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Number of clusters") +
    ggplot2::ylab("Within-cluster sum of distances")
}
