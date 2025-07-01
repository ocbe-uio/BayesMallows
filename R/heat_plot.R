#' Heat plot of posterior probabilities
#'
#' Generates a heat plot with items in their consensus ordering along the
#' horizontal axis and ranking along the vertical axis. The color denotes
#' posterior probability.
#'
#' @param model_fit An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param ... Additional arguments passed on to other methods. In particular,
#'   `type = "CP"` or `type = "MAP"` can be passed on to
#'   [compute_consensus()] to determine the order of items along the
#'   horizontal axis.
#'
#' @return A ggplot object.
#' @export
#' @details In models with a single cluster, the items are sorted along the
#'   x-axis according to the CP consensus ranking. In models with more than one
#'   cluster, the items are sorted alphabetically.
#'
#' @example /inst/examples/heat_plot_example.R
#' @family posterior quantities
heat_plot <- function(model_fit, ...) {
  if (is.null(burnin(model_fit))) {
    stop("Please specify the burnin with 'burnin(model_fit) <- value'.")
  }

  if(model_fit$n_clusters == 1) {
    item_order <- unique(compute_consensus(model_fit, ...)[["item"]])
  } else {
    item_order <- model_fit$data$items
  }


  posterior_ranks <- model_fit$rho[
    model_fit$rho$iteration > burnin(model_fit), ,
    drop = FALSE
  ]
  posterior_ranks$probability <- 1

  heatplot_data <- aggregate(posterior_ranks[, "probability", drop = FALSE],
    by = list(
      cluster = posterior_ranks$cluster,
      item = posterior_ranks$item,
      value = posterior_ranks$value
    ),
    FUN = function(x) sum(x) / length(unique(posterior_ranks$iteration))
  )

  heatplot_data$item <- factor(heatplot_data$item, levels = item_order)
  heatplot_data <- heatplot_data[order(heatplot_data$item), , drop = FALSE]

  heatplot_expanded <- expand.grid(
    cluster = unique(heatplot_data$cluster),
    item = unique(heatplot_data$item),
    value = unique(heatplot_data$value)
  )
  heatplot_expanded <- merge(
    heatplot_expanded, heatplot_data,
    by = c("cluster", "item", "value"), all.x = TRUE
  )
  heatplot_expanded$probability[is.na(heatplot_expanded$probability)] <- 0

  p <- ggplot2::ggplot(
    heatplot_expanded,
    ggplot2::aes(x = .data$item, y = .data$value, fill = .data$probability)
  ) +
    ggplot2::geom_tile() +
    ggplot2::labs(fill = "Probability") +
    ggplot2::xlab("Item") +
    ggplot2::ylab("Rank")

  if(model_fit$n_clusters > 1) {
    p + ggplot2::facet_wrap(ggplot2::vars(.data$cluster))
  } else {
    p
  }
}
