#' Posterior plots
#'
#' @param fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_posterior}}.
#' @param type Which type of plot. Either \code{"alpha"} for a posterior histogram
#'   of the scale parameter or \code{"rho"} for a posterior histogram of a ranked
#'   item.
#' @param item If \code{type = "rho"} this argument must be set. It defines
#'   which item to plot. If \code{type = "alpha"} this argument is ignored.
#' @param ... Other arguments passed to ggplot2::geom_histogram.
#'
#' @return A plot.
#' @export
#'
#' @examples
#' # Fit a model to the potato data
#' model_fit <- compute_posterior(potato_weighing, "footrule",
#' nmc = 10000, burnin = 5000)
#' # Plot the posterior histogram of alpha
#' plot(model_fit, type = "alpha", bins = 50)
plot.BayesMallows <- function(fit, type = "alpha", item = NULL, ...){
  if(type == "alpha") {
    df <- data.frame(alpha = fit$alpha)

    ggplot2::ggplot(df, ggplot2::aes(x = alpha)) +
      ggplot2::geom_histogram(...) +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ggtitle(paste("Mean = ", round(mean(df$alpha), 2),
                             "standard deviation =", round(sd(df$alpha), 2)))




  }
}
