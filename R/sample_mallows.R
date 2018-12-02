#' Random Samples from the Mallows Rank Model
#'
#' Generate random samples from the Mallows Rank Model
#' \insertCite{mallows1957}{BayesMallows} with consensus ranking \eqn{\rho} and
#' scale parameter \eqn{\alpha}. The samples are obtained by running the
#' Metropolis-Hastings algorithm described in Appendix C of
#' \insertCite{vitelli2018;textual}{BayesMallows}.
#'
#' @param rho0 Vector specifying the latent consensus ranking in the Mallows
#'   rank model.
#' @param alpha0 Scalar specifying the scale parameter in the Mallows rank
#'   model.
#' @param n_samples Integer specifying the number of random samples to generate.
#'   When \code{diagnostic = TRUE}, this number must be larger than 1.
#' @param leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution.
#' @param metric Character string specifying the distance measure to use.
#'   Available options are \code{"footrule"} (default), \code{"spearman"},
#'   \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and \code{"ulam"}.
#'   See also the \code{rmm} function in the \code{PerMallows} package
#'   \insertCite{irurozki2016}{BayesMallows} for sampling from the Mallows
#'   model with Cayley, Hamming, Kendall, and Ulam distances.
#' @param diagnostic Logical specifying whether to output convergence
#'   diagnostics. If \code{TRUE}, a diagnostic plot is printed, together with
#'   the returned samples.
#' @param burnin Integer specifying the number of iterations to discard as
#'   burn-in. Defaults to 1000 when \code{diagnostic = FALSE}, else to 0.
#' @param thinning Integer specifying the number of MCMC iterations to perform
#'   between each time a random rank vector is sampled. Defaults to 1000 when
#'   \code{diagnostic = FALSE}, else to 1.
#' @param items_to_plot Integer vector used if \code{diagnostic = TRUE}, in
#'   order to specify the items to plot in the diagnostic output. If not
#'   provided, 5 items are picked at random.
#' @param max_lag Integer specifying the maximum lag to use in the computation
#'   of autocorrelation. Defaults to 1000L. This argument is passed to
#'   \code{stats::acf}. Only used when \code{diagnostic = TRUE}.
#'
#' @references \insertAllCited{}
#'
#' @export
#'
#' @example /inst/examples/sample_mallows_example.R
#'
sample_mallows <- function(rho0, alpha0, n_samples,
                           leap_size = 1,
                           metric = "footrule",
                           diagnostic = FALSE,
                           burnin = ifelse(diagnostic, 0, 1000),
                           thinning = ifelse(diagnostic, 1, 1000),
                           items_to_plot = NULL,
                           max_lag = 1000L)
                            {

  if(!(validate_permutation(rho0) && sum(is.na(rho0)) == 0)){
    stop("rho0 must be a proper ranking with no missing values.")
  }

  if(diagnostic && n_samples == 1){
    stop("Must have more than one samples to create diagnostic plots")
  } else if(n_samples <= 0){
    stop("n_samples must be positive.")
  }

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
    if(is.null(items_to_plot) && n_items > 5){
      message("Items not provided by user. Picking 5 at random.")
      items_to_plot <- sample.int(n_items, 5)
    } else {
      items_to_plot <- seq(from = 1, to = n_items, by = 1)
    }

    # Compute the autocorrelation in the samples
    autocorr <- apply(samples[ , items_to_plot, drop = FALSE], 2, stats::acf,
                      lag.max = max_lag, plot = FALSE, demean = TRUE)
    names(autocorr) <- items_to_plot


    autocorr <- purrr::map_dfr(autocorr, function(x) {
      dplyr::tibble(
        acf = x$acf[, 1, 1],
        lag = x$lag[, 1, 1]
      )
    }, .id = "item")
    autocorr <- dplyr::mutate(autocorr, item = as.factor(as.integer(.data$item)))

    ac_plot <- ggplot2::ggplot(autocorr,
                               ggplot2::aes(x = .data$lag, y = .data$acf, color = .data$item)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Lag") +
      ggplot2::ylab("Autocorrelation") +
      ggplot2::ggtitle("Autocorrelation of Rank Values")

    colnames(samples) <- seq(from = 1, to = n_items, by = 1)
    diagnostic <- dplyr::as_tibble(samples)
    diagnostic <- dplyr::mutate(diagnostic, iteration = dplyr::row_number())

    diagnostic <- tidyr::gather(diagnostic, key = "item",
                                value = "value", -.data$iteration)
    diagnostic <- dplyr::filter(diagnostic, .data$item %in% items_to_plot)
    diagnostic <- dplyr::mutate(diagnostic, item = as.factor(as.integer(.data$item)))



    rho_plot <- ggplot2::ggplot(diagnostic,
                         ggplot2::aes(x = .data$iteration, y = .data$value, color = .data$item)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.title = ggplot2::element_blank()) +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab("Rank value") +
      ggplot2::ggtitle("Trace Plot of Rank Values")

    print(cowplot::plot_grid(ac_plot, rho_plot, ncol = 1))

    samples <- samples[seq(from = burnin + 1, by = thinning, length.out = n_samples), ]
  }
  return(samples)

}
