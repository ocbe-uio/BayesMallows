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
#' options are \code{"alpha"} and \code{"rho"}.
#'
#' @param items Numeric vector specifying the index of the items to plot
#' or character vector with names of the items, as specified in
#'   \code{model_fit$items}. Only used when \code{type = "rho"}.
#'
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @export
#'
#' @examples
#' # Analysis of complete rankings
#' # The example datasets potato_visual and potato_weighing contain complete
#' # rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
#' # model:
#' model_fit <- compute_mallows(potato_visual)
#' # We study the trace plot of the parameters
#' # alpha is the default
#' assess_convergence(model_fit)
#' # When studying convergence of rho, we can also specify which items to plot
#' assess_convergence(model_fit, type = "rho", items = 1:5)
#' # Based on these plots, we conclude that the Markov chain has converged well
#' # before 1,000 iterations. We hence set burnin = 1000.
#' # Next, we use the generic plot function to study the posterior distributions
#' # of alpha and rho
#' plot(model_fit, burnin = 1000)
#' plot(model_fit, burnin = 1000, type = "rho", items = 1:20)
#'
plot.BayesMallows <- function(x, burnin, type = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.

  stopifnot(x$nmc > burnin)

  stopifnot(type %in% c("alpha", "rho"))

  if(type == "alpha") {
    df <- dplyr::filter(x$alpha, .data$iteration > burnin)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_density() +
      ggplot2::xlab(expression(alpha)) +
      ggplot2::ylab("Posterior density") +
      ggplot2::ggtitle(label = "Posterior density of alpha")

    if(x$n_clusters > 1){
      p <- p + ggplot2::facet_wrap(~ .data$cluster)
    }

    print(p)

  } else if(type == "rho") {

    if(is.null(items)) stop("You must specify the items to plot.")

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
  }
}
