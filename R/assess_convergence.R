#' Trace Plots from Metropolis-Hastings Algorithm
#'
#' \code{assess_convergence} provides trace plots for the parameters of the
#' Mallows Rank model, in order to study the convergence of the Metropolis-Hastings
#' algorithm.
#'
#' @param model_fit A fitted model object of class \code{BayesMallows}, obtained
#'   with \code{\link{compute_mallows}}.
#'
#' @param type Character string specifying which plot type we want. Available
#' options are \code{"alpha"}, \code{"rho"}, \code{"Rtilde"}, or
#' \code{"cluster_probs"}.
#'
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{model_fit$items} or a
#'   vector of indices. If NULL, five items are selected randomly. Only used when \code{type = "rho"}.
#'
#' @param assessors Numeric vector specifying the assessors to study in
#' the diagnostic plot for \code{"Rtilde"}.
#'
#' @seealso \code{\link{compute_mallows}}, \code{\link{plot.BayesMallows}}
#'
#' @export
#'
assess_convergence <- function(model_fit, type = "alpha", items = NULL,
                               assessors = NULL){

  stopifnot(class(model_fit) == "BayesMallows")

  if(type == "alpha") {

    # Create the diagnostic plot for alpha
    p <- ggplot2::ggplot(model_fit$alpha, ggplot2::aes(x = .data$iteration, y = .data$value)) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(expression(alpha)) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::ggtitle(label = "Convergence of alpha")

    if(model_fit$n_clusters == 1){
      p <- p + ggplot2::geom_line()
    } else {
      p <- p + ggplot2::geom_line(ggplot2::aes(color = .data$cluster))
    }

    print(p)

  } else if(type == "rho"){

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
      ggplot2::ylab(expression(rho)) +
      ggplot2::ggtitle(label = "Convergence of rho")

    if(model_fit$n_clusters > 1){
      p <- p + ggplot2::facet_wrap(~ .data$cluster)
    }

    print(p)

  } else if(type == "Rtilde") {

    if(!model_fit$save_augmented_data){
      stop("Please rerun with compute_mallows with save_augmented_data = TRUE")
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
    }

    if(is.factor(model_fit$augmented_data$item) && is.numeric(items)){
      items <- levels(model_fit$augmented_data$item)[items]
    }
    df <- dplyr::filter(model_fit$augmented_data,
                        .data$assessor %in% assessors,
                        .data$item %in% items)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~ .data$assessor) +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Rtilde") +
      ggplot2::ggtitle(label = "Convergence of Rtilde")


  } else if (type == "cluster_probs"){

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
    stop("type must be either \"alpha\", \"rho\", \"augmentation\", or \"cluster_probs\".")
  }
}
