#' Generate elbow plot for choice of clusters
#'
#' @param ... One or more objects returned from \code{\link{compute_mallows}}.
#' @param rankings A matrix with complete rankings, the same as was
#' given to \code{\link{compute_mallows}}.
#' @param burnin The number of iterations to discard.
#' @param probs Probability of the quantiles to returned. This argument is passed on to
#' \code{stats::quantile}.
#' @param thinning Optional argument which can be used to take only every
#' \code{thinning}th iteration after burn-in. If the function takes very long to run,
#' one might try to set this number to something larger than 1. Be aware that this always leads to
#' less accurate posterior distributions.
#'
#' @details At the moment, this function computes the footrule distance,
#' and no other distance measures are supported.
#'
#' @export
#'
plot_elbow <- function(..., rankings, burnin,
                       probs = seq(0, 1, 0.25),
                       thinning = 1){

  # Put the models into a list. These are typically fitted with different number of clusters
  models <- list(...)

  # Compute within-cluster distance for each, and put in a dataframe
  distances <- purrr::map_dfr(models, .compute_within_cluster_distance,
                              rankings = rankings, burnin = burnin, thinning = thinning)


  # Do a boxplot
  ggplot2::ggplot(distances, ggplot2::aes(x = as.factor(.data$n_clusters), y = .data$within_cluster_distance)) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Number of clusters") +
    ggplot2::ylab("Within-cluster sum of distances") +
    ggplot2::ggtitle("Elbow Plot")

}
