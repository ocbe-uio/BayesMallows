#
#' Compute Posterior Intervals
#'
#' Compute posterior intervals of parameters of interest.
#'
#' @param model_fit A model object.
#' @param parameter Character string defining which parameter to compute
#'   posterior intervals for. One of `"alpha"`, `"rho"`, or
#'   `"cluster_probs"`. Default is `"alpha"`.
#' @param level Decimal number in \eqn{[0,1]} specifying the confidence level.
#'   Defaults to `0.95`.
#' @param decimals Integer specifying the number of decimals to include in
#'   posterior intervals and the mean and median. Defaults to `3`.
#' @param ... Other arguments. Currently not used.
#'
#' @details This function computes both the Highest Posterior Density Interval (HPDI),
#' which may be discontinuous for bimodal distributions, and
#' the central posterior interval, which is simply defined by the quantiles of the posterior
#' distribution.
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/compute_posterior_intervals_example.R
#'
#' @export
#' @family posterior quantities
compute_posterior_intervals <- function(model_fit, ...) {
  UseMethod("compute_posterior_intervals")
}

#' @export
#' @rdname compute_posterior_intervals
compute_posterior_intervals.BayesMallows <- function(
    model_fit, parameter = c("alpha", "rho", "cluster_probs"),
    level = 0.95, decimals = 3L, ...) {
  if (is.null(burnin(model_fit))) {
    stop("Please specify the burnin with 'burnin(model_fit) <- value'.")
  }

  parameter <- match.arg(parameter, c("alpha", "rho", "cluster_probs"))

  stopifnot(level > 0 && level < 1)

  posterior_data <- model_fit[[parameter]][
    model_fit[[parameter]]$iteration > burnin(model_fit), ,
    drop = FALSE
  ]

  if (parameter == "alpha" || parameter == "cluster_probs") {
    posterior_split <- split(posterior_data, f = posterior_data$cluster)

    posterior_intervals <- do.call(rbind, lapply(posterior_split, function(x) {
      data.frame(
        parameter = parameter,
        cluster = unique(x$cluster),
        mean = format(round(mean(x$value), decimals), nsmall = decimals),
        median = format(round(stats::median(x$value), decimals),
          nsmall = decimals
        ),
        hpdi = compute_continuous_hpdi(x$value, level, decimals),
        central_interval = compute_central_interval(x$value, level, decimals)
      )
    }))
  } else if (parameter == "rho") {
    posterior_split <- split(
      posterior_data,
      f = list(posterior_data$item, posterior_data$cluster)
    )

    posterior_intervals <- do.call(rbind, lapply(posterior_split, function(x) {
      data.frame(
        parameter = parameter,
        cluster = unique(x$cluster),
        item = unique(x$item),
        mean = round(mean(x$value), 0),
        median = round(stats::median(x$value), 0),
        hpdi = compute_discrete_hpdi(x, level),
        central_interval = compute_central_interval(x$value, level, 0)
      )
    }))
  }

  if (model_fit$n_clusters == 1) posterior_intervals$cluster <- NULL
  row.names(posterior_intervals) <- NULL
  posterior_intervals
}

#' @export
#' @rdname compute_posterior_intervals
compute_posterior_intervals.SMCMallows <- function(
    model_fit, parameter = c("alpha", "rho"), level = 0.95,
    decimals = 3L, ...) {
  model_fit$compute_options$burnin <- 0
  parameter <- match.arg(parameter, c("alpha", "rho"))
  NextMethod("compute_posterior_intervals")
}

compute_central_interval <- function(values, level, decimals) {
  central <- unique(
    stats::quantile(values,
      probs = c((1 - level) / 2, level + (1 - level) / 2),
      names = FALSE
    )
  )
  interval <- format(round(central, decimals), nsmall = decimals)
  paste0("[", paste(trimws(interval), collapse = ","), "]")
}

# This function is derived from HDInterval::hdiVector
# Copyright: Juat Ngumbang, Mike Meredith, and John Kruschke
compute_continuous_hpdi <- function(values, level, decimals) {
  n <- length(values)
  values <- sort(values)
  lower <- values[1:(n - floor(n * level))]
  upper <- values[(floor(n * level) + 1):n]
  ind <- which.min(upper - lower)
  hpdi <- format(round(c(lower[ind], upper[ind]), decimals), nsmall = decimals)
  paste0("[", paste(trimws(hpdi), collapse = ","), "]")
}

compute_discrete_hpdi <- function(x, level) {
  pct_dat <- aggregate(
    iteration ~ value,
    data = x, FUN = function(y) {
      length(y) / nrow(x)
    }
  )
  pct_dat <- pct_dat[order(pct_dat$iteration, decreasing = TRUE), ]
  pct_dat$cumprob <- cumsum(pct_dat$iteration)
  maxind <- min(which(pct_dat$cumprob >= level))
  hpdi <- sort(pct_dat$value[seq(from = 1, to = maxind)])
  contiguous_regions <- split(hpdi, cummax(c(1, diff(hpdi))))
  hpdi <- vapply(contiguous_regions, function(r) {
    paste0("[", paste(unique(range(r)), collapse = ","), "]")
  }, character(1))
  paste(hpdi, collapse = "")
}
