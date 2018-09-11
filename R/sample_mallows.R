#' Random Samples from the Mallows Rank Model
#'
#' Generate random samples from the Mallows
#' Rank Model \insertCite{mallows1957}{BayesMallows} with consensus ranking \eqn{\rho}
#' and scale parameter \eqn{\alpha}. The samples are obtained by running the
#' Metropolis-Hastings algorithm described in Appendix C of
#' \insertCite{vitelli2018;textual}{BayesMallows}.
#'
#' @param rho0 Vector specifying the latent consensus ranking.
#' @param alpha0 Scalar specifying the scale parameter.
#' @param n_samples Integer specyfing the number of random samples to generate.
#' @param burnin Integer specifying the number of iterations to discard as burn-in.
#' @param thinning Integer specifying the number of MCMC iterations to perform
#' between each time a random rank vector is sampled.
#' @param leap_size Integer specifying the step size of the leap-and-shift proposal distribution.
#' @param metric Character string specifying the distance measure to use. Available
#' options are \code{"footrule"} (default), \code{"spearman"}, \code{"cayley"}, and
#' \code{"kendall"}. For sampling from the Mallows model with Cayley and Kendall distances
#' the \code{PerMallows} package \insertCite{irurozki2016}{BayesMallows} can also be used.
#' @param diagnostic Logical specifying whether to output convergence diagnostics. If \code{TRUE},
#' a diagnostic plot is returned together with the samples.
#' @param items_to_plot Integer vector used if \code{diagnostic = TRUE}, in order to
#' specify the items to plot in the diagnostic output. If not provided, 5 items are picked
#' at random.
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @example /inst/examples/sample_mallows_example.R
#'
sample_mallows <- function(rho0, alpha0, n_samples,
                           burnin, thinning, leap_size = 1,
                           metric = "footrule",
                           diagnostic = FALSE,
                           items_to_plot = NULL)
                            {
  n_items <- length(rho0)

  if(diagnostic){
    internal_burnin <- 0
    internal_thinning <- 1
    internal_n_samples <- burnin + n_samples * thinning
  } else {
    internal_burnin <- burnin
    internal_thinning <- thinning
    internal_n_samples <- n_samples
  }

  samples <- t(rmallows(
    rho0 = rho0,
    alpha0 = alpha0,
    n_samples = internal_n_samples,
    burnin = internal_burnin,
    thinning = internal_thinning,
    leap_size = leap_size,
    metric = metric
  ))

  if(diagnostic){
    if(is.null(items_to_plot)){
      message("Items not provided by user. Picking 5 at random.")
      items_to_plot <- sample.int(n_items, 5)
    }

    diagnostic <- dplyr::as_tibble(samples)
    names(diagnostic) <- seq(from = 1, to = n_items, by = 1)
    diagnostic <- dplyr::mutate(diagnostic, iteration = dplyr::row_number())

    diagnostic <- tidyr::gather(diagnostic, key = "item", value = "value", -.data$iteration)
    diagnostic <- dplyr::filter(diagnostic, .data$item %in% items_to_plot)
    diagnostic <- dplyr::mutate(diagnostic,
                                item = as.factor(as.integer(.data$item)))

    p <- ggplot2::ggplot(diagnostic,
            ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Sampled value")

    print(p)

    samples <- samples[seq(from = burnin + 1, by = thinning, length.out = n_samples), ]
  }
  return(samples)

}
