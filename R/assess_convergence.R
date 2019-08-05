#' Trace Plots from Metropolis-Hastings Algorithm
#'
#' \code{assess_convergence} provides trace plots for the parameters of the
#' Mallows Rank model, in order to study the convergence of the Metropolis-Hastings
#' algorithm.
#'
#' @param model_fit A fitted model object of class \code{BayesMallows} returned from
#'  \code{\link{compute_mallows}} or an object of class \code{BayesMallowsMixtures}
#'  returned from \code{\link{compute_mallows_mixtures}}.
#'
#' @param parameter Character string specifying which parameter to plot. Available
#' options are \code{"alpha"}, \code{"rho"}, \code{"Rtilde"},
#' \code{"cluster_probs"}, or \code{"theta"}.
#'
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{model_fit$items} or a
#'   vector of indices. If NULL, five items are selected randomly. Only used when
#'   \code{parameter = "rho"} or \code{parameter = "Rtilde"}.
#'
#' @param assessors Numeric vector specifying the assessors to study in
#' the diagnostic plot for \code{"Rtilde"}.
#'
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @export
#'
assess_convergence <- function(model_fit, parameter = "alpha", items = NULL,
                               assessors = NULL){

  stopifnot(inherits(model_fit, "BayesMallows") ||
              inherits(model_fit, "BayesMallowsMixtures"))

  if(parameter == "alpha") {
    if(inherits(model_fit, "BayesMallows")){
      m <- model_fit$alpha
      trace_alpha(m, FALSE)
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      m <- purrr::map_dfr(model_fit, function(x){
        dplyr::mutate(x$alpha,
                      cluster = as.character(.data$cluster),
                      n_clusters = x$n_clusters)
      })
      trace_alpha(m, TRUE)
      }

  } else if(parameter == "rho"){
    if(inherits(model_fit, "BayesMallows")){
      trace_rho(model_fit, items)
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      cowplot::plot_grid(plotlist = purrr::map(model_fit, trace_rho, clusters = TRUE, items = items))
    }

  } else if(parameter == "Rtilde") {

    if(inherits(model_fit, "BayesMallows")){
      trace_rtilde(model_fit, items, assessors)
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      stop("Trace plots of augmented data not supported for BayesMallowsMixtures. Please rerun each component k using the k-th list element.")
    }
  } else if (parameter == "cluster_probs"){
    if(inherits(model_fit, "BayesMallows")){
      m <- model_fit$cluster_probs
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      m <- purrr::map_dfr(model_fit, function(x){
        dplyr::mutate(x$cluster_probs,
                      cluster = as.character(.data$cluster),
                      n_clusters = x$n_clusters)
      })
    }
    trace_cluster_probs(m)

  } else if (parameter == "theta"){
      trace_theta(model_fit)
  } else {
    stop("parameter must be either \"alpha\", \"rho\", \"augmentation\", \"cluster_probs\", or \"theta\".")
  }
}

trace_alpha <- function(m, clusters){
  # Create the diagnostic plot for alpha
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data$iteration, y = .data$value)) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha))

  if(!clusters){
    p <- p + ggplot2::geom_line()
  } else {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(color = .data$cluster)) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::facet_wrap(ggplot2::vars(.data$n_clusters),
                          labeller = ggplot2::as_labeller(function(n_clusters) {
                            paste(n_clusters, dplyr::if_else(n_clusters == 1, "cluster", "clusters"))
                          }), scales = "free_y")
  }
  return(p)
}

trace_rho <- function(model_fit, items, clusters = model_fit$n_clusters > 1){

  if(is.null(items) && model_fit$n_items > 5){
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(model_fit$n_items, 5)
  } else if (is.null(items) && model_fit$n_items > 0) {
    items <- seq.int(from = 1, to = model_fit$n_items)
  }

  if(!is.character(items)){
    items <- model_fit$items[items]
  }

  df <- dplyr::filter(model_fit$rho, .data$item %in% items)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(rho))

  if(clusters){
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$cluster))
  }

  return(p)
}


trace_rtilde <- function(model_fit, items, assessors, ...){


  if(!model_fit$save_aug){
    stop("Please rerun with compute_mallows with save_aug = TRUE")
  }

  if(is.null(items) && model_fit$n_items > 5){
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(model_fit$n_items, 5)
  } else if (is.null(items) && model_fit$n_items > 0) {
    items <- seq.int(from = 1, to = model_fit$n_items)
  }

  if(is.null(assessors) && model_fit$n_assessors > 5){
    message("Assessors not provided by user. Picking 5 at random.")
    assessors <- sample.int(model_fit$n_assessors, 5)
  } else if (is.null(assessors) && model_fit$n_assessors > 0) {
    assessors <- seq.int(from = 1, to = model_fit$n_assessors)
  } else if(!is.null(assessors)) {
    if(length(setdiff(assessors, seq(1, model_fit$n_assessors, 1))) > 0) {
      stop("assessors vector must contain numeric indices between 1 and the number of assessors")
    }
  }

  if(is.factor(model_fit$augmented_data$item) && is.numeric(items)){
    items <- levels(model_fit$augmented_data$item)[items]
  }
  df <- dplyr::filter(model_fit$augmented_data,
                      .data$assessor %in% assessors,
                      .data$item %in% items)
  df <- dplyr::mutate(df, assessor = as.factor(.data$assessor))
  levels(df$assessor) <- paste("Assessor", levels(df$assessor))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(ggplot2::vars(.data$assessor)) +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab("Rtilde")
}


trace_cluster_probs <- function(m){

  ggplot2::ggplot(m, ggplot2::aes(x = .data$iteration, y = .data$value,
                               color = .data$cluster)) +
    ggplot2::geom_line() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(tau[c])) +
    ggplot2::facet_wrap(ggplot2::vars(.data$n_clusters),
                        labeller = ggplot2::as_labeller(function(n_clusters) {
                          paste(n_clusters, dplyr::if_else(n_clusters == 1, "cluster", "clusters"))
                        }), scales = "free_y")
}


trace_theta <- function(model_fit){
  if(is.null(model_fit$theta)){
    stop("Theta not available. Run compute_mallows with error_model = 'bernoulli'.")
  }
  # Create the diagnostic plot for theta
  p <- ggplot2::ggplot(model_fit$theta, ggplot2::aes(x = .data$iteration, y = .data$value)) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(theta)) +
    ggplot2::geom_line()

  return(p)
}
