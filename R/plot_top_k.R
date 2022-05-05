#' Plot Top-k Rankings with Pairwise Preferences
#'
#' Plot the posterior probability, per item, of being ranked among the top-\eqn{k}
#' for each assessor. This plot is useful when the data take the form of pairwise
#' preferences.
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#'
#' @param k Integer specifying the k in top-\eqn{k}.
#'
#' @param rel_widths The relative widths of the plots of \code{rho} per cluster
#' and the plot of assessors, respectively. This argument is passed on to
#' \code{\link[cowplot]{plot_grid}}.
#'
#' @seealso \code{\link{predict_top_k}}
#'
#' @export
#'
#' @example /inst/examples/plot_top_k_example.R
#'
plot_top_k <- function(model_fit, burnin = model_fit$burnin,
                       k = 3,
                       rel_widths = c(model_fit$n_clusters, 10)) {
  validate_top_k(model_fit, burnin)

  # Extract post burn-in rows with value <= k
  rho <- model_fit$rho[model_fit$rho$iteration > burnin & model_fit$rho$value <= k, , drop = FALSE]
  n_samples <- length(unique(rho$iteration))
  # Factors are not needed in this case
  rho$item <- as.character(rho$item)
  rho <- aggregate(
    list(prob = rho$iteration),
    list(item = rho$item, cluster = rho$cluster),
    FUN = function(x) length(x) / n_samples
  )

  # Find the complete set of items per cluster
  rho <- do.call(rbind, lapply(split(rho, f = rho$cluster), function(dd) {
    dd2 <- merge(dd, expand.grid(item = unique(rho$item)),
      by = "item", all = TRUE
    )
    dd2$cluster[is.na(dd2$cluster)] <- unique(dd$cluster)
    dd2$prob[is.na(dd2$prob)] <- 0
    dd2
  }))[, c("cluster", "item", "prob")]

  # Sort the items according to probability in Cluster 1
  item_ordering <- compute_consensus(model_fit, type = "CP", burnin = burnin)
  if (model_fit$n_clusters > 1) {
    item_ordering <- rev(item_ordering[item_ordering$cluster == "Cluster 1", ]$item)
  } else {
    item_ordering <- rev(item_ordering$item)
  }

  rho$item <- factor(rho$item, levels = unique(item_ordering))

  # Trick to make the plot look nicer
  if (model_fit$n_clusters == 1) {
    rho$cluster <- ""
  }

  rankings <- .predict_top_k(model_fit, burnin = burnin, k = k)

  # Sorting the items according to their probability in rho
  rankings$item <- factor(rankings$item, levels = item_ordering)

  assessor_plot <- ggplot2::ggplot(rankings, ggplot2::aes(.data$assessor, .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "white") +
    ggplot2::xlab("Assessor") +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )

  rho_plot <- ggplot2::ggplot(rho, ggplot2::aes(.data$cluster, .data$item)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$prob), colour = "white") +
    ggplot2::ylab("Item") +
    ggplot2::xlab(expression(rho)) +
    ggplot2::theme(legend.position = "none")

  if (model_fit$n_clusters > 1) {
    rho_plot <- rho_plot + ggplot2::facet_wrap(~ .data$cluster)
  }

  cowplot::plot_grid(rho_plot, assessor_plot, rel_widths = rel_widths)
}
