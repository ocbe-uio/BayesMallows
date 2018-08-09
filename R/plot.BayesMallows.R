#' @title Posterior plots
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#' @param burnin Number of iterations to discard as burn-in.
#' @param type Either \code{"alpha"} or \code{"rho"}.
#' @param items Index of items to plot. Only used when \code{type = "rho"}.
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @return A plot.
#' @export
#'
plot.BayesMallows <- function(x, burnin, type = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.

  stopifnot(x$nmc > burnin)

  stopifnot(type %in% c("alpha", "rho"))


  if(type == "alpha") {
    start <- floor(burnin / x$alpha_jump) + 1
    stop <- ceiling(x$nmc / x$alpha_jump)
    num <- stop - start + 1

    df <- dplyr::tibble(alpha = x$alpha[start:stop, 1])

    ggplot2::ggplot(df, ggplot2::aes(x = .data$alpha)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density") +
      ggplot2::ggtitle(label = "Posterior density of alpha",
                       subtitle = paste(num, "samples, mean =",
                       sprintf("%.1f", mean(df$alpha))))

  } else if(type == "rho") {
    if(is.null(items)) stop("You must specify the items to plot.")

    start <- floor(burnin / x$thinning) + 1
    stop <- ceiling(x$nmc / x$thinning)

    df <- gather_rho(x, items, row_inds = start:stop)

    # Compute the density, rather than the count, since the latter
    # depends on the number of Monte Carlo samples
    df <- dplyr::group_by(df, .data$Item, .data$Rank)
    df <- dplyr::summarise(df, n = dplyr::n())
    df <- dplyr::mutate(df, pct = .data$n / sum(.data$n))

    # Function for getting an x axis without decimals.
    # Modified from https://stackoverflow.com/questions/21061653/creating-a-density-histogram-in-ggplot2
    scalefun <- function(x) sprintf("%d", as.integer(x))

    # Finally create the plot
    ggplot2::ggplot(df, ggplot2::aes(x = .data$Rank, y = .data$pct)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(labels = scalefun) +
      ggplot2::facet_wrap(~ .data$Item) +
      ggplot2::ggtitle("Posterior ranks for items") +
      ggplot2::xlab("Rank") +
      ggplot2::ylab("Posterior probability")
  }
}
