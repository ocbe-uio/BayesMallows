#' @title Posterior plots
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @return A plot.
#' @export
#'
plot.BayesMallows <- function(x, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.

  graphics::plot(x$alpha, type = "l")
  # if(type == "alpha") {
  #   df <- dplyr::tibble(alpha = model_fit$alpha)
  #
  #   ggplot2::ggplot(df, ggplot2::aes_(x =~ alpha)) +
  #     ggplot2::geom_histogram(...) +
  #     ggplot2::xlab(expression(alpha)) +
  #     ggplot2::ggtitle(paste("Mean = ", round(mean(df$alpha), 2),
  #                            "standard deviation =", round(stats::sd(df$alpha), 2)))
  #
  # } else if(type == "rho") {
  #   if(is.null(items)) stop("You must specify the items to plot.")
  #
  #   df <- gather_rho(model_fit, items)
  #
  #   # Compute the density, rather than the count, since the latter
  #   # depends on the number of Monte Carlo samples
  #   df <- dplyr::group_by(df, .data$Item, .data$Rank)
  #   df <- dplyr::summarise(df, n = n())
  #   df <- dplyr::mutate(df, pct = n / sum(n))
  #
  #   # Function for getting an x axis without decimals.
  #   # Modified from https://stackoverflow.com/questions/21061653/creating-a-density-histogram-in-ggplot2
  #   scalefun <- function(x) sprintf("%d", as.integer(x))
  #
  #   # Finally create the plot
  #   ggplot2::ggplot(df, ggplot2::aes(x = .data$Rank, y = .data$pct)) +
  #     ggplot2::geom_col(...) +
  #     ggplot2::scale_x_continuous(labels = scalefun) +
  #     ggplot2::facet_wrap(~ .data$Item) +
  #     ggplot2::ggtitle("Posterior ranks for items") +
  #     ggplot2::xlab("Rank") +
  #     ggplot2::ylab("Posterior probability")
  #
  #
  # }
}
