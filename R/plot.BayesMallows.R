#' Posterior plots
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#' @param type Which type of plot. Either \code{"alpha"} for a posterior histogram
#'   of the scale parameter or \code{"rho"} for a posterior histogram of a ranked
#'   item.
#' @param items If \code{type = "rho"} this argument must be set. It defines
#'   which item(s) to plot. If \code{type = "alpha"} this argument is ignored.
#' @param ... Other arguments passed to ggplot2::geom_histogram.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' # Fit a model to the potato data
#' model_fit <- compute_mallows(potato_weighing, "footrule",
#' nmc = 10000, burnin = 5000)
#' # Plot the posterior histogram of alpha
#' plot(model_fit, type = "alpha", bins = 50)
#' # We can also plot the posterior histograms of rankings of items
#' plot(model_fit, type = "rho", items = c(2, 5), bins = 20)
#'
plot.BayesMallows <- function(x, type = "alpha", items = NULL, ...){
  if(type == "alpha") {
    df <- dplyr::tibble(alpha = x$alpha)

    ggplot2::ggplot(df, ggplot2::aes_(x =~ alpha)) +
      ggplot2::geom_histogram(...) +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ggtitle(paste("Mean = ", round(mean(df$alpha), 2),
                             "standard deviation =", round(stats::sd(df$alpha), 2)))

  } else if(type == "rho") {
    if(is.null(items)) stop("You must specify the items to plot.")

    # Convert the posterior rank matrix to a tibble
    df <- dplyr::as_tibble(t(x$rho))

    # Save the number of items
    n <- ncol(df)

    # Set the column names to the item numbers
    names(df) <- seq(1, n)

    # Make the tibble tall by gathering items
    df <- tidyr::gather(df, key = "Item", value = "Rank")

    # Filter the tibble to keep the items the user wants
    df <- dplyr::filter(df, .data$Item %in% items)

    # Compute the density, rather than the count, since the latter
    # depends on the number of Monte Carlo samples
    df <- dplyr::group_by(df, .data$Item, .data$Rank)
    df <- dplyr::summarise(df, n = n())
    df <- dplyr::mutate(df, pct = n / sum(n))

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
