#
#' Compute Posterior Intervals
#'
#' Compute posterior intervals of parameters of interest.
#'
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#'
#' @param parameter Character string defining which parameter to compute
#' posterior intervals for. One of \code{"alpha"}, \code{"rho"}, or
#' \code{"cluster_probs"}. Default is \code{"alpha"}.
#'
#' @param level Decimal number in \eqn{[0,1]} specifying the confidence level.
#' Defaults to \code{0.95}.
#'
#' @param decimals Integer specifying the number of decimals to include
#' in posterior intervals and the mean and median. Defaults to \code{3}.
#'
#' @details This function computes both the Highest Posterior Density Interval (HPDI),
#' which may be discontinuous for bimodal distributions, and
#' the central posterior interval, which is simply defined by the quantiles of the posterior
#' distribution. The HPDI intervals are computed using the \code{HDInterval} package
#' \insertCite{meredith2018}{BayesMallows}.
#'
#' @seealso \code{\link{compute_mallows}}
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/compute_posterior_intervals_example.R
#'
#' @export
#'
compute_posterior_intervals <- function(model_fit, burnin = model_fit$burnin,
                                        parameter = "alpha", level = 0.95,
                                        decimals = 3L){
  stopifnot(class(model_fit) == "BayesMallows")

  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)
  stopifnot(parameter %in% c("alpha", "rho", "cluster_probs", "cluster_assignment"))
  stopifnot(level > 0 && level < 1)

  df <- dplyr::filter(model_fit[[parameter]], .data$iteration > burnin)

  if(parameter == "alpha" || parameter == "cluster_probs"){

    df <- dplyr::group_by(df, .data$cluster)
    df <- .compute_posterior_intervals(df, parameter, level, decimals)

  } else if(parameter == "rho"){
    decimals <- 0
    df <- dplyr::group_by(df, .data$cluster, .data$item)
    df <- .compute_posterior_intervals(df, parameter, level, decimals, discrete = TRUE)

  }

  df <- dplyr::ungroup(df)

  if(model_fit$n_clusters == 1) df <- dplyr::select(df, -.data$cluster)

  return(df)
}


.compute_posterior_intervals <- function(df, parameter, level, decimals, discrete = FALSE){
  dplyr::do(df, {
    format <- paste0("%.", decimals, "f")

    posterior_mean <- round(base::mean(.data$value), decimals)
    posterior_median <- round(stats::median(.data$value), decimals)

    if(discrete) {

      df <- dplyr::group_by(.data, .data$value)
      df <- dplyr::summarise(df, n = dplyr::n())
      df <- dplyr::arrange(df, dplyr::desc(.data$n))
      df <- dplyr::mutate(df, cumprob = cumsum(.data$n) / sum(.data$n),
                          lagcumprob = dplyr::lag(.data$cumprob, default = 0))

      df <- dplyr::filter(df, .data$lagcumprob <= level)

      values <- sort(dplyr::pull(df, .data$value))

      # Find contiguous regions
      breaks <- c(0, which(diff(values) != 1), length(values))

      hpdi <- purrr::map(seq(length(breaks) - 1), function(.x, values, breaks) {
        vals <- values[(breaks[.x] + 1):breaks[.x + 1]]
        vals <- unique(c(min(vals), max(vals)))
        paste0("[", paste(vals, collapse = ","), "]")
        }, values = values, breaks = breaks)

      hpdi <- paste(hpdi, collapse = ",")

    } else {
      hpdi <- HDInterval::hdi(.data$value, credMass = level, allowSplit = TRUE)

      hpdi[] <- sprintf(format, hpdi)
      if(is.matrix(hpdi)){
        # Discontinous case
        hpdi <- paste(apply(hpdi, 1, function(x) paste0("[", x[[1]], ",", x[[2]], "]")))
      } else {
        # Continuous case
        hpdi <- paste0("[", hpdi[[1]], ",", hpdi[[2]], "]")
      }
    }


    central <- unique(stats::quantile(.data$value, probs = c((1 - level) / 2, level + (1 - level) / 2)))
    central <- sprintf(format, central)
    central <- paste0("[", paste(central, collapse = ","), "]")

    dplyr::tibble(
      parameter = parameter,
      mean = posterior_mean,
      median = posterior_median,
      conf_level = paste(level * 100, "%"),
      hpdi = hpdi,
      central_interval = central
    )
  })
}
