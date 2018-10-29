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
#' options are \code{"alpha"}, \code{"rho"}, \code{"Rtilde"}, or
#' \code{"cluster_probs"}.
#'
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{model_fit$items} or a
#'   vector of indices. If NULL, five items are selected randomly. Only used when
#'   \code{parameter = "rho"} or \code{parameter = "Rtilde"}.
#'
#' @param assessors Numeric vector specifying the assessors to study in
#' the diagnostic plot for \code{"Rtilde"}.
#'
#' @param ... Additional arguments passed on to \code{cowplot::plot_grid}. Only used
#' when \code{model_fit} is of class \code{BayesMallowsMixtures}.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @export
#'
assess_convergence <- function(model_fit, parameter = "alpha", items = NULL,
                               assessors = NULL, ...){

  stopifnot(inherits(model_fit, "BayesMallows") ||
              inherits(model_fit, "BayesMallowsMixtures"))


  if(parameter == "alpha") {
    if(inherits(model_fit, "BayesMallows")){
      trace_alpha(model_fit)
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      cowplot::plot_grid(plotlist = purrr::map(model_fit, trace_alpha, clusters = TRUE), ...)
    }

  } else if(parameter == "rho"){
    if(inherits(model_fit, "BayesMallows")){
      trace_rho(model_fit, items)
    } else if(inherits(model_fit, "BayesMallowsMixtures")){
      cowplot::plot_grid(plotlist = purrr::map(model_fit, trace_rho, clusters = TRUE, items = items))
    }



  } else if(parameter == "Rtilde") {

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
      ggplot2::facet_wrap(~ .data$assessor) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Rtilde") +
      ggplot2::ggtitle(label = "Convergence of Rtilde")


  } else if (parameter == "cluster_probs"){

    if(!exists("cluster_probs", model_fit)){
      stop("cluster_probs not found")
    }

    ggplot2::ggplot(model_fit$cluster_probs,
                    ggplot2::aes(x = .data$iteration, y = .data$value,
                                 color = .data$cluster)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(tau[k])) +
      ggplot2::ggtitle("Cluster Probabilities")

  } else {
    stop("parameter must be either \"alpha\", \"rho\", \"augmentation\", or \"cluster_probs\".")
  }
}



trace_alpha <- function(model_fit, clusters = model_fit$n_clusters > 1){
  # Create the diagnostic plot for alpha
  p <- ggplot2::ggplot(model_fit$alpha, ggplot2::aes(x = .data$iteration, y = .data$value)) +
    ggplot2::xlab("Iteration") +
    ggplot2::ylab(expression(alpha))

  if(!clusters){
    p <- p + ggplot2::geom_line()
  } else {
    p <- p +
      ggplot2::geom_line(ggplot2::aes(color = .data$cluster)) +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::theme(legend.title = ggplot2::element_blank())
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
    p <- p + ggplot2::facet_wrap(~ .data$cluster)
  }

  return(p)
}
