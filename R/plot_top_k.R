
#' Find the top-k probabilities
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#' @param burnin Number of iterations to discard as burn-in.
#' @param k The k in top-k
#' @param rel_widths The relative widths of the plots of \code{rho} per cluster
#' and the plot of assessors, respectively. This argument is passed on to
#' \code{cowplot::plot_grid}.
#'
#' @export
#'
plot_top_k <- function(model_fit, burnin, k = 3, rel_widths = c(rep(1, model_fit$n_clusters), 10)){
  stopifnot(burnin < model_fit$nmc)

  n_samples <- sum(unique(model_fit$rho$iteration) > burnin)

  # Extract post burn-in rows with value <= k
  rho <- dplyr::filter(model_fit$rho, .data$iteration > burnin, .data$value <= k)
  # Factors are not needed in this case
  rho <- dplyr::mutate(rho, item = as.character(.data$item))
  rho <- dplyr::group_by(rho, .data$item, .data$cluster)
  rho <- dplyr::summarise(rho, prob = dplyr::n()/n_samples)
  rho <- dplyr::ungroup(rho)

  # Find the complete set of items per cluster
  rho <- tidyr::complete(
    dplyr::group_by(rho, .data$cluster),
    item = model_fit$items,
    fill = list(prob = 0))
  rho <- dplyr::ungroup(rho)

  # Sort the items according to probability
  item_ordering <- rev(compute_cp_consensus(model_fit, burnin = burnin)$item)
  rho <- dplyr::mutate(rho, item = factor(.data$item, levels = item_ordering))

  # Trick to make the plot look nicer
  if(model_fit$n_clusters == 1){
    rho <- dplyr::mutate(rho, cluster = "")
  }

  rankings <- dplyr::filter(model_fit$augmented_data, .data$iteration > burnin, .data$value <= k)
  rankings <- dplyr::mutate(rankings, item = as.character(.data$item))
  rankings <- dplyr::group_by(rankings, .data$assessor, .data$item)
  rankings <- dplyr::summarise(rankings, prob = dplyr::n()/n_samples)
  rankings <- dplyr::ungroup(rankings)
  rankings <- tidyr::complete(
    dplyr::group_by(rankings, .data$assessor),
    item = model_fit$items,
    fill = list(prob = 0)
  )
  # Sorting the items according to their probability in rho
  rankings <- dplyr::mutate(rankings, item = factor(.data$item, levels = item_ordering))

  assessor_plot <- ggplot2::ggplot(rankings, ggplot2::aes(.data$assessor, .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "white") +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::xlab("Assessor") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
      )


  rho_plot <- ggplot2::ggplot(rho, ggplot2::aes(.data$cluster, .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "white") +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::ylab("Item") +
    ggplot2::xlab(expression(rho)) +
    ggplot2::theme(legend.position = "none")

  if(model_fit$n_clusters > 1){
    rho_plot <- rho_plot + ggplot2::facet_wrap(~ .data$cluster)
  }

  cowplot::plot_grid(rho_plot, assessor_plot, rel_widths = rel_widths)
}
