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
#' @param burnin The number of iterations to discard as burnin. The following
#' options are possible:
#' \itemize{
#' \item A numeric vector with one element for each model object provided in
#' \code{...}.
#'
#' \item A single number. In this case, the burnin for each model provided in
#' \code{...} is taken to be equal to this single number.
#'
#' \item If the argument to \code{...} is a list of models, and this list
#' has an element named \code{burnin}, then this is the default value of \code{burnin}.
#' The element can either be a vector or a single number, as described in the points
#' above.
#'
#' \item If the arguments to \code{...} are a set of models, and each
#' models has a \code{burnin} element set, then this is the default value of
#' \code{burnin}.
#' }
#'
#' @return A boxplot with the number of clusters on the horizontal axis and the
#' with-cluster sum of distances on the vertical axis.
#'
#' @export
#'
#' @seealso \code{\link{compute_mallows}}
#'
plot_elbow <- function(..., burnin = NULL){

  # Put the models into a list. These are typically fitted with different number of clusters
  models <- list(...)

  # Taking into account the case where the user has entered a list of models
  if(length(models) == 1 && !(class(models[[1]]) == "BayesMallows")){
    models <- models[[1]]
  }


  #### TODO: Deal with the case where burnin is either
  # an element of models or an element of the elements inside models.

  # Set the burnin
  if(!is.null(burnin)){
    if(length(burnin) == 1 && length(models) > 1){
      models <- purrr::map(models, function(x) {
        x$burnin <- burnin
        return(x)
      })
    } else if(length(burnin) == length(models)){
      models <- purrr::map2(models, burnin, function(x, y){
        x$burnin <- y
        return(x)
      })
    } else {
      stop("Burnin must be either a single number or a vector of same length
           as the number of models.")
    }
  }

  df <- purrr::map_dfr(models, function(x) {
    stopifnot(class(x) == "BayesMallows")
    if(!x$include_wcd) stop("To get an elbow plot, set include_wcd=TRUE in compute_mallows")

    df <- dplyr::filter(x$within_cluster_distance, .data$iteration > x$burnin)

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
