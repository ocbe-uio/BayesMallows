#' @title Posterior plots
#'
#' @param x An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#' @param burnin Number of iterations to discard as burn-in.
#' @param type Either \code{"alpha"} or \code{"rho"}.
#' @param items Index of items to plot or names of the items, as specified in
#'   \code{rownames(model_fit$rho)}. Only used when \code{type = "rho"}.
#' @param ... Other arguments passed to \code{plot} (not used).
#'
#' @return A plot.
#' @export
#'
#' @examples
#' # Here is an example.
plot.BayesMallows <- function(x, burnin, type = "alpha", items = NULL, ...){
  # Note, the first argument must be named x, otherwise R CMD CHECK will
  # issue a warning. This is because plot.BayesMallows must have the same
  # required arguments as graphics::plot.

  stopifnot(x$nmc > burnin)

  stopifnot(type %in% c("alpha", "rho"))



  if(type == "alpha") {
    start <- floor(burnin / x$alpha_jump) + 1
    stop <- ceiling(x$nmc / x$alpha_jump)

    alpha_matrix <- x$alpha[seq(from = start, to = stop, by = 1), , drop = FALSE]

    df <- prepare_alpha_df(alpha_matrix)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$alpha)) +
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

    start <- floor(burnin / x$thinning) + 1
    stop <- ceiling(x$nmc / x$thinning)

    df <- gather_rho(x, items, row_inds = seq(from = start, to = stop, by = 1))

    # Compute the density, rather than the count, since the latter
    # depends on the number of Monte Carlo samples
    df <- dplyr::group_by(df, .data$cluster, .data$item, .data$rank)
    df <- dplyr::summarise(df, n = dplyr::n())
    df <- dplyr::mutate(df, pct = .data$n / sum(.data$n))

    # Finally create the plot
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$rank, y = .data$pct)) +
      ggplot2::geom_col() +
      ggplot2::scale_x_continuous(labels = scalefun) +
      #ggplot2::scale_y_continuous(limits = c(0, 1)) +
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
