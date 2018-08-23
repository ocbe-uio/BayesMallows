#' @title Posterior plots
#'
#' @param ... One or more objects returned from \code{\link{compute_mallows}},
#' separated by comma, or a list of such objects.
#' @param burnin Number of iterations to discard as burn-in.
#'
#' @return A plot.
#' @export
#'
plot_elbow <- function(..., burnin){

  # Put the models into a list. These are typically fitted with different number of clusters
  models <- list(...)

  # Taking into account the case where the user has entered a list of models
  if(length(models) == 1 && !(class(models[[1]]) == "BayesMallows")){
    models <- models[[1]]
  }

  df <- purrr::map_dfr(models, function(x) {
    stopifnot(class(x) == "BayesMallows")
    if(!x$include_wcd) stop("To get an elbow plot, set include_wcd=TRUE in compute_mallows")

    df <- dplyr::filter(x$within_cluster_distance, .data$iteration > burnin)

    # Need to sum the within-cluster distances across clusters, for each iteration
    df <- dplyr::group_by(df, .data$iteration)
    df <- dplyr::summarise(df, value = sum(.data$value))

    df <- dplyr::mutate(df, n_clusters = x$n_clusters)
    return(df)
  })

  ggplot2::ggplot(df, ggplot2::aes(x = as.factor(.data$n_clusters), y = .data$value)) +
    ggplot2::geom_boxplot() +
    ggplot2::xlab("Number of clusters") +
    ggplot2::ylab("Within-cluster som of distances") +
    ggplot2::ggtitle("Elbow Plot")


}
