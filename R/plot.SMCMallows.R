#' @title Plot SMC Posterior Distributions
#' @description Plot posterior distributions of SMC-Mallow parameters.
#' @param x An object of type `SMC-Mallows`, returned for example from
#' [smc_mallows_new_users()].
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to `model_fit$burnin`, and must be
#' provided if `model_fit$burnin` does not exist. See
#' [assess_convergence()].
#' @param parameter Character string defining the parameter to plot. Available
#' options are `"alpha"` and `"rho"`.
#' @param time Integer determining the update slice to plot
#' @param C Number of cluster
#' @param colnames A vector of item names. If NULL, generic names are generated
#' for the items in the ranking.
#' @param items Either a vector of item names, or a vector of indices. If NULL,
#' five items are selected randomly.
#' @param ... Other arguments passed to [base::plot()] (not used).
#' @return A plot of the posterior distributions
#' @author Waldir Leoncio
#' @export
#' @example /inst/examples/plot.SMCMallows_example.R
#' @family posterior quantities
plot.SMCMallows <- function(
    x, nmc = nrow(x$rho_samples[, 1, ]), burnin = 0,
    parameter = "alpha", time = ncol(x$rho_samples[, 1, ]), C = 1,
    colnames = NULL, items = NULL, ...) {
  if (parameter == "alpha") {
    output <- x$alpha_samples[, time]
    plot_alpha_smc(output, nmc, burnin)
  } else if (parameter == "rho") {
    output <- x$rho_samples[, , time]
    plot_rho_smc(output, nmc, burnin, C, colnames, items)
  } else {
    stop("parameter must be either 'alpha' or 'rho'.")
  }
}

plot_alpha_smc <- function(output, nmc, burnin) {
  alpha_samples_table <- data.frame(iteration = 1:nmc, value = output)

  ggplot2::ggplot(alpha_samples_table, ggplot2::aes(x = .data$value)) +
    ggplot2::geom_density() +
    ggplot2::xlab(expression(alpha)) +
    ggplot2::ylab("Posterior density") +
    ggplot2::ggtitle(label = "Implemented SMC scheme") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

plot_rho_smc <- function(output, nmc, burnin, C, colnames = NULL, items = NULL) {
  n_items <- dim(output)[2]

  if (is.null(items) && n_items > 5) {
    message("Items not provided by user or more than 5 items in a ranking. Picking 5 at random.")
    items <- sample(seq_len(n_items), 5, replace = FALSE)
    items <- sort(items)
  } else if (is.null(items) && n_items <= 5) {
    items <- seq_len(n_items)
    items <- sort(items)
  }

  # do smc processing here
  smc_plot <- smc_processing(output = output, colnames = colnames)

  if (!is.character(items)) {
    items <- unique(smc_plot$item)[items]
  }

  iteration <- rep(seq_len(nmc), times = n_items)
  df <- cbind(iteration, smc_plot)

  if (C == 1) {
    df <- cbind(cluster = "Cluster 1", df)
  }

  df <- df[df$iteration > burnin & df$item %in% items, , drop = FALSE]

  # Compute the density, rather than the count, since the latter
  # depends on the number of Monte Carlo samples
  df <- aggregate(list(n = df$iteration),
    list(cluster = df$cluster, item = df$item, value = df$value),
    FUN = length
  )
  df$pct <- df$n / sum(df$n)

  df$item <- factor(df$item, levels = c(items))

  # Finally create the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pct)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(labels = scalefun) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("Posterior probability")

  if (C == 1) {
    p <- p + ggplot2::facet_wrap(~ .data$item)
  } else {
    p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
  }

  return(p)
}
