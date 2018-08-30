#' Plot Posterior Distributions
#'
#' Plot posterior distributions of the parameters of the Mallows Rank model.
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. See \code{\link{assess_convergence}}.
#'
#' @param type Character string defining the parameter to plot. Available
#' options are \code{"alpha"}, \code{"rho"}, \code{"cluster_probs"}, and
#' \code{"cluster_assignment"}.
#'
#' @param items The items to study in the diagnostic plot for \code{rho}. Either
#'   a vector of item names, corresponding to \code{x$items} or a
#'   vector of indices. If NULL, five items are selected randomly.
#'   Only used when \code{type = "rho"}.
#'
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @export
#'
plot.BayesMallows <- function(x, burnin, type = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.

  stopifnot(x$nmc > burnin)

  stopifnot(type %in% c("alpha", "rho", "cluster_probs", "cluster_assignment"))

  if(type == "alpha") {
    df <- dplyr::filter(x$alpha, .data$iteration > burnin)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density") +
      ggplot2::ggtitle(label = "Posterior density of alpha")

    if(x$n_clusters > 1){
      p <- p + ggplot2::facet_wrap(~ .data$cluster, scales = "free_x")
    }

    print(p)

  } else if(type == "rho") {

    if(is.null(items) && x$n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items <- sample.int(x$n_items, 5)
    } else if (is.null(items) && x$n_items > 0) {
      items <- seq.int(from = 1, to = x$n_items)
    }

    if(!is.character(items)){
      items <- x$items[items]
    }

    df <- dplyr::filter(x$rho, .data$iteration > burnin, .data$item %in% items)

    # Compute the density, rather than the count, since the latter
    # depends on the number of Monte Carlo samples
    df <- dplyr::group_by(df, .data$cluster, .data$item, .data$value)
    df <- dplyr::summarise(df, n = dplyr::n())
    df <- dplyr::mutate(df, pct = .data$n / sum(.data$n))

    # Finally create the plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pct)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(labels = scalefun) +
      ggplot2::ggtitle("Posterior ranks for items") +
      ggplot2::xlab("rank") +
      ggplot2::ylab("Posterior probability")

    if(x$n_clusters == 1){
      p <- p + ggplot2::facet_wrap(~ .data$item)
    } else {
      p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
    }

    print(p)
  } else if(type == "cluster_probs"){
    df <- dplyr::filter(x$cluster_probs, .data$iteration > burnin)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(tau[k])) +
      ggplot2::ylab("Posterior density") +
      ggplot2::ggtitle(label = "Posterior density of cluster probabilities") +
      ggplot2::facet_wrap(~ .data$cluster)

  } else if(type == "cluster_assignment"){

    # First get one cluster per assessor, and sort these
    df <- assign_cluster(x, burnin = burnin, soft = FALSE, expand = FALSE)
    df <- dplyr::arrange(df, .data$map_cluster)
    assessor_order <- dplyr::pull(df, .data$assessor)

    # Next, rerun with soft=TRUE to get the probability of all clusters
    df <- assign_cluster(x, burnin = burnin, soft = TRUE, expand = TRUE)
    # Then order the assessors according to assessor_order
    df <- dplyr::mutate(df, assessor = factor(.data$assessor, levels = assessor_order))

    # Now make a plot
    ggplot2::ggplot(df, ggplot2::aes(.data$assessor, .data$cluster)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$probability)) +
      ggplot2::scale_fill_gradient(low = "blue", high = "red") +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      ) +
      ggplot2::xlab(paste0("Assessors (", min(assessor_order), " - ", max(assessor_order), ")")) +
      ggplot2::ggtitle("Posterior Probabilities of Cluster Assignment")

  }
}
