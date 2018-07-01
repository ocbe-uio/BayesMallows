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
    df <- data.frame(alpha = x$alpha)

    ggplot2::ggplot(df, ggplot2::aes_(x =~ alpha)) +
      ggplot2::geom_histogram(...) +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ggtitle(paste("Mean = ", round(mean(df$alpha), 2),
                             "standard deviation =", round(stats::sd(df$alpha), 2)))

  } else if(type == "rho") {
    stopifnot(!is.null(items))

    nmc <- length(x$rho_acceptance)
    df <- data.frame(
      Item = as.factor(items),
      Rank = as.factor(x$rho[items, ])
      )

    ggplot2::ggplot(df, ggplot2::aes_(x =~ Rank)) +
      ggplot2::geom_histogram(stat = "count", ...) +
      ggplot2::facet_wrap(~ Item) +
      ggplot2::ggtitle(paste("Posterior histograms of items", paste(items, collapse = ", ")))
  }
}
