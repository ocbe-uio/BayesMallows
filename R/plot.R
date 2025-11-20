#' Plot Posterior Distributions
#'
#' Plot posterior distributions of the parameters of the Mallows Rank model.
#'
#' @param x An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param parameter Character string defining the parameter to plot. Available
#' options are `"alpha"`, `"rho"`, `"cluster_probs"`,
#' `"cluster_assignment"`, and `"theta"`.
#'
#' @param items The items to study in the diagnostic plot for `rho`. Either
#'   a vector of item names, corresponding to `x$data$items` or a
#'   vector of indices. If NULL, five items are selected randomly.
#'   Only used when `parameter = "rho"`.
#'
#' @param ... Other arguments passed to `plot` (not used).
#'
#' @export
#' @importFrom stats density
#'
#' @example /inst/examples/plot.BayesMallows_example.R
#' @family posterior quantities
plot.BayesMallows <- function(x, parameter = "alpha", items = NULL, ...) {
  parameter <- match.arg(
    parameter,
    c("alpha", "rho", "cluster_probs", "cluster_assignment", "theta")
  )

  if (is.null(burnin(x))) {
    stop("Please specify the burnin with 'burnin(x) <- value'.")
  }

  if (parameter == "alpha") {
    plot_alpha(x)
  } else if (parameter == "rho") {
    plot_rho(x, items)
  } else if (parameter == "cluster_probs") {
    df <- x$cluster_probs[x$cluster_probs$iteration > burnin(x), , drop = FALSE]

    ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(tau[c])) +
      ggplot2::ylab("Posterior density") +
      ggplot2::facet_wrap(~ .data$cluster)
  } else if (parameter == "cluster_assignment") {
    if (is.null(x$cluster_assignment)) {
      stop("No cluster assignments.")
    }

    df <- assign_cluster(x, soft = FALSE, expand = FALSE)
    df <- df[order(df$map_cluster), ]
    assessor_order <- df$assessor

    df <- assign_cluster(x, soft = TRUE, expand = TRUE)
    df$assessor <- factor(df$assessor, levels = assessor_order)

    ggplot2::ggplot(df, ggplot2::aes(.data$assessor, .data$cluster)) +
      ggplot2::geom_tile(ggplot2::aes(fill = .data$probability)) +
      ggplot2::theme(
        legend.title = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank()
      ) +
      ggplot2::xlab(paste0("Assessors (", min(assessor_order), " - ", max(assessor_order), ")"))
  } else if (parameter == "theta") {
    if (is.null(x$theta)) {
      stop("Please run compute_mallows with error_model = 'bernoulli'.")
    }
    df <- x$theta[x$theta$iteration > burnin(x), , drop = FALSE]
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(theta)) +
      ggplot2::ylab("Posterior density")
    return(p)
  }
}


#' @title Plot SMC Posterior Distributions
#' @description Plot posterior distributions of SMC-Mallow parameters.
#' @param x An object of type \code{SMC-Mallows}.
#' @param parameter Character string defining the parameter to plot. Available
#'   options are \code{"alpha"} and \code{"rho"}.
#' @param items Either a vector of item names, or a vector of indices. If NULL,
#'   five items are selected randomly.
#' @param ... Other arguments passed to \code{\link[base]{plot}} (not used).
#' @return A plot of the posterior distributions
#' @export
#' @family posterior quantities
#' @example /inst/examples/update_mallows_example.R
#'
plot.SMCMallows <- function(
    x, parameter = "alpha", items = NULL, ...) {
  parameter <- match.arg(parameter, c("alpha", "rho"))

  if (parameter == "alpha") {
    plot_alpha(x)
  } else if (parameter == "rho") {
    plot_rho(x, items)
  } else {
    stop("parameter must be either 'alpha' or 'rho'.")
  }
}

plot_alpha <- function(x) {
  plot_dat <- x$alpha[x$alpha$iteration > burnin(x), , drop = FALSE]

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density() +
    ggplot2::xlab(expression(alpha)) +
    ggplot2::ylab("Posterior density")

  if (x$n_clusters > 1) {
    p <- p + ggplot2::facet_wrap(~ .data$cluster, scales = "free_x")
  }
  p
}


plot_rho <- function(x, items) {
  if (is.null(items) && x$data$n_items > 5) {
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(x$data$n_items, 5)
  } else if (is.null(items) && x$data$n_items > 0) {
    items <- seq.int(from = 1, to = x$data$n_items)
  } else if (!is.null(items)) {
    if (!all(items %in% x$data$items) && !all(items %in% seq_along(x$data$items))) {
      stop("Unknown items.")
    }
  }

  if (!is.character(items)) {
    items <- x$data$items[items]
  }

  df <- x$rho[x$rho$iteration > burnin(x) & x$rho$item %in% items, , drop = FALSE]
  df1 <- aggregate(iteration ~ item + cluster + value, data = df, FUN = length)
  df1$pct <- df1$iteration / length(unique(df$iteration)) / length(unique(df$chain))

  # Finally create the plot
  p <- ggplot2::ggplot(df1, ggplot2::aes(x = factor(.data$value), y = .data$pct)) +
    ggplot2::geom_col() +
    ggplot2::xlab("rank") +
    ggplot2::ylab("Posterior probability")

  if (x$n_clusters == 1) {
    p <- p + ggplot2::facet_wrap(~ .data$item)
  } else {
    p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
  }

  p
}
