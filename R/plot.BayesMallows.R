#' Plot Posterior Distributions
#'
#' Plot posterior distributions of the parameters of the Mallows Rank model.
#'
#' @param x An object of type `BayesMallows`, returned from
#'   [compute_mallows()].
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to `x$burnin`, and must be
#' provided if `x$burnin` does not exist. See [assess_convergence()].
#'
#' @param parameter Character string defining the parameter to plot. Available
#' options are `"alpha"`, `"rho"`, `"cluster_probs"`,
#' `"cluster_assignment"`, and `"theta"`.
#'
#' @param items The items to study in the diagnostic plot for `rho`. Either
#'   a vector of item names, corresponding to `x$items` or a
#'   vector of indices. If NULL, five items are selected randomly.
#'   Only used when `parameter = "rho"`.
#'
#' @param ... Other arguments passed to `plot` (not used).
#'
#' @export
#'
#' @example /inst/examples/plot.BayesMallows_example.R
#' @family posterior quantities
plot.BayesMallows <- function(x, burnin = x$burnin, parameter = "alpha", items = NULL, ...) {

  parameter <- match.arg(
    parameter,
    c("alpha", "rho", "cluster_probs", "cluster_assignment", "theta"))

  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  if (x$nmc <= burnin) stop("burnin must be <= nmc")

  if (parameter == "alpha") {

    plot_alpha(x, burnin)

  } else if (parameter == "rho") {
    plot_rho(x, items, burnin)
  } else if (parameter == "cluster_probs") {
    df <- x$cluster_probs[x$cluster_probs$iteration > burnin, , drop = FALSE]

    ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(tau[c])) +
      ggplot2::ylab("Posterior density") +
      ggplot2::facet_wrap(~ .data$cluster)
  } else if (parameter == "cluster_assignment") {
    if (is.null(x$cluster_assignment)) {
      stop("No cluster assignments.")
    }

    # First get one cluster per assessor, and sort these
    df <- assign_cluster(x, burnin = burnin, soft = FALSE, expand = FALSE)
    df <- df[order(df$map_cluster), ]
    assessor_order <- df$assessor

    # Next, rerun with soft=TRUE to get the probability of all clusters
    df <- assign_cluster(x, burnin = burnin, soft = TRUE, expand = TRUE)
    # Then order the assessors according to assessor_order
    df$assessor <- factor(df$assessor, levels = assessor_order)

    # Now make a plot
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

    df <- x$theta[x$theta$iteration > burnin, , drop = FALSE]

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

plot_alpha <- function(x, burnin = 0) {
  plot_dat <- x$alpha[x$alpha$iteration > burnin, , drop = FALSE]

  p <- ggplot2::ggplot(plot_dat, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density() +
    ggplot2::xlab(expression(alpha)) +
    ggplot2::ylab("Posterior density")

  if (x$n_clusters > 1) {
    p <- p + ggplot2::facet_wrap(~ .data$cluster, scales = "free_x")
  }
  p
}


plot_rho <- function(x, items, burnin = 0) {
  if (is.null(items) && x$n_items > 5) {
    message("Items not provided by user. Picking 5 at random.")
    items <- sample.int(x$n_items, 5)
  } else if (is.null(items) && x$n_items > 0) {
    items <- seq.int(from = 1, to = x$n_items)
  }

  if (!is.character(items)) {
    items <- x$items[items]
  }

  df <- x$rho[x$rho$iteration > burnin & x$rho$item %in% items, , drop = FALSE]

  # Compute the density, rather than the count, since the latter
  # depends on the number of Monte Carlo samples
  df <- aggregate(
    list(n = df$iteration),
    list(cluster = df$cluster, item = df$item, value = df$value),
    FUN = length
  )
  df$pct <- df$n / sum(df$n)

  # Finally create the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pct)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(labels = scalefun) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("Posterior probability")

  if (x$n_clusters == 1) {
    p <- p + ggplot2::facet_wrap(~ .data$item)
  } else {
    p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
  }

  return(p)
}
