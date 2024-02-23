#' Plot Within-Cluster Sum of Distances
#'
#' Plot the within-cluster sum of distances from the corresponding cluster
#' consensus for different number of clusters. This function is useful for
#' selecting the number of mixture.
#'
#' @param ... One or more objects returned from [compute_mallows()],
#' separated by comma, or a list of such objects. Typically, each object
#' has been run with a different number of mixtures, as specified in the
#' `n_clusters` argument to [compute_mallows()].
#'
#' @return A boxplot with the number of clusters on the horizontal axis and the
#' with-cluster sum of distances on the vertical axis.
#'
#' @export
#'
#' @example /inst/examples/compute_mallows_mixtures_example.R
#' @family posterior quantities
plot_elbow <- function(...) {
  # Put the models into a list. These are typically fitted with different number of clusters
  models <- list(...)

  # Taking into account the case where the user has entered a list of models
  if (length(models) == 1 && !(inherits(models[[1]], "BayesMallows"))) {
    models <- models[[1]]
  }

  df <- do.call(rbind, lapply(models, function(x) {
    stopifnot(inherits(x, "BayesMallows"))

    if (is.null(burnin(x))) {
      stop("Please specify burnin with 'burnin(model_fit) <- value'.")
    }

    if (length(unique(x$within_cluster_distance$iteration)) != x$compute_options$nmc) {
      stop("To get an elbow plot, set include_wcd=TRUE in compute_mallows")
    }

    df <- x$within_cluster_distance[
      x$within_cluster_distance$iteration > burnin(x), ,
      drop = FALSE
    ]

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
