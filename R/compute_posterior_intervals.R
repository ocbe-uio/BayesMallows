#
#' Compute Posterior Intervals
#'
#' Compute posterior intervals of parameters of interest.
#'
#' @param model_fit A model object.
#'
#' @param ... other arguments passed to methods.
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
#' @family posterior quantities
compute_posterior_intervals <- function(model_fit, ...) {
  UseMethod("compute_posterior_intervals")
}


#' Compute posterior intervals
#'
#' @param model_fit An object of class \code{BayesMallows} returned from
#'   \code{\link{compute_mallows}}.
#' @param burnin A numeric value specifying the number of iterations to discard
#'   as burn-in. Defaults to \code{model_fit$burnin}, and must be provided if
#'   \code{model_fit$burnin} does not exist. See
#'   \code{\link{assess_convergence}}.
#' @param parameter Character string defining which parameter to compute
#'   posterior intervals for. One of \code{"alpha"}, \code{"rho"}, or
#'   \code{"cluster_probs"}. Default is \code{"alpha"}.
#' @param level Decimal number in \eqn{[0,1]} specifying the confidence level.
#'   Defaults to \code{0.95}.
#' @param decimals Integer specifying the number of decimals to include in
#'   posterior intervals and the mean and median. Defaults to \code{3}.
#' @param ... Other arguments. Currently not used.
#'
#' @seealso assess_convergence
#' @export
#' @family posterior quantities
compute_posterior_intervals.BayesMallows <- function(
    model_fit, burnin = model_fit$burnin,
    parameter = "alpha",
    level = 0.95, decimals = 3L, ...) {
  stopifnot(inherits(model_fit, "BayesMallows"))

  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)

  if (length(parameter) > 1) stop("Only one parameter allowed.")
  parameter <- match.arg(
    parameter,
    c("alpha", "rho", "cluster_probs", "cluster_assignment")
  )

  stopifnot(level > 0 && level < 1)

  df <- model_fit[[parameter]][model_fit[[parameter]]$iteration > burnin, , drop = FALSE]
  df$chain <- NULL

  if (parameter == "alpha" || parameter == "cluster_probs") {
    df <- .compute_posterior_intervals(split(df, f = df$cluster), parameter, level, decimals)
  } else if (parameter == "rho") {
    decimals <- 0
    df <- .compute_posterior_intervals(
      split(df, f = interaction(df$cluster, df$item)),
      parameter, level, decimals,
      discrete = TRUE
    )
  }

  if (model_fit$n_clusters == 1) df$cluster <- NULL

  row.names(df) <- NULL
  return(df)
}

#' @title Compute posterior intervals
#'
#' @description This function computes posterior intervals based on the set of samples at the
#' last timepoint of the SMC algorithm.
#'
#' @param model_fit An object of class \code{SMCMallows}, returned from
#'   \code{\link{smc_mallows_new_item_rank}} or
#'   \code{\link{smc_mallows_new_users}}.
#' @param parameter Character string defining which parameter to compute
#'   posterior intervals for. One of \code{"alpha"} or \code{"rho"}.
#' @param level Decimal number in \eqn{[0,1]} specifying the confidence level.
#'   Defaults to \code{0.95}.
#' @param decimals Integer specifying the number of decimals to include in
#'   posterior intervals and the mean and median. Defaults to \code{3}.
#' @param ... Other arguments. Currently not used.
#' @export
#' @family posterior quantities
#'
#' @example inst/examples/smc_post_processing_functions_example.R
compute_posterior_intervals.SMCMallows <- function(
    model_fit, parameter = "alpha", level = 0.95,
    decimals = 3L, ...) {
  if (length(parameter) > 1) stop("Only one parameter allowed.")
  parameter <- match.arg(
    parameter, c("alpha", "rho")
  )

  stopifnot(level > 0 && level < 1)

  if (parameter == "alpha") {
    tab <- data.frame(
      value = model_fit$alpha_samples[, ncol(model_fit$alpha_samples), drop = TRUE]
    )
    tab$n_clusters <- 1
    tab$cluster <- "Cluster 1"

    tab <- .compute_posterior_intervals(
      df = split(tab, f = tab$cluster),
      parameter = parameter,
      level = level,
      decimals = decimals
    )
  } else if (parameter == "rho") {
    tab <- smc_processing(
      model_fit$rho_samples[, , dim(model_fit$rho_samples)[[3]], drop = TRUE]
    )
    tab$n_clusters <- 1
    tab$cluster <- "Cluster 1"

    tab <- .compute_posterior_intervals(
      df = split(tab, f = interaction(tab$cluster, tab$item)),
      parameter = parameter, level = level,
      decimals = decimals, discrete = TRUE
    )
  }

  if (length(unique(tab$cluster)) == 1) {
    tab$cluster <- NULL
  }

  rownames(tab) <- NULL

  return(tab)
}

.compute_posterior_intervals <- function(df, parameter, level, decimals, discrete = FALSE, ...) {
  do.call(rbind, lapply(df, function(x) {
    format <- paste0("%.", decimals, "f")

    posterior_mean <- round(base::mean(x$value), decimals)
    posterior_median <- round(stats::median(x$value), decimals)

    if (discrete) {
      hpdi <- compute_discrete_hpdi(x, level)
    } else {
      hpdi <- HDInterval::hdi(x$value, credMass = level, allowSplit = TRUE)

      hpdi[] <- sprintf(format, hpdi)
      if (is.matrix(hpdi)) {
        # Discontinous case
        hpdi <- paste(apply(hpdi, 1, function(x) paste0("[", x[[1]], ",", x[[2]], "]")))
      } else {
        # Continuous case
        hpdi <- paste0("[", hpdi[[1]], ",", hpdi[[2]], "]")
      }
    }

    central <- unique(stats::quantile(x$value, probs = c((1 - level) / 2, level + (1 - level) / 2)))
    central <- sprintf(format, central)
    central <- paste0("[", paste(central, collapse = ","), "]")

    ret <- data.frame(
      parameter = parameter,
      mean = posterior_mean,
      median = posterior_median,
      conf_level = paste(level * 100, "%"),
      hpdi = hpdi,
      central_interval = central
    )

    targets <- setdiff(names(x), c("iteration", "value", "n_clusters"))
    for (nm in targets) {
      eval(parse(text = paste0("ret$", nm, " <- unique(x$", nm, ")")))
    }
    ret[, c(targets, setdiff(names(ret), targets)), drop = FALSE]
  }))
}


compute_discrete_hpdi <- function(df, level) {
  if (!"iteration" %in% names(df)) df$iteration <- seq_len(nrow(df))
  df <- aggregate(list(n = df$iteration),
    list(value = df$value),
    FUN = length
  )
  df <- df[order(df$n, decreasing = TRUE), , drop = FALSE]
  df$cumprob <- cumsum(df$n) / sum(df$n)
  df$lagcumprob <- c(0, head(df$cumprob, -1))
  df <- df[df$lagcumprob <= level, , drop = FALSE]

  values <- sort(df$value)

  # Find contiguous regions
  breaks <- c(0, which(diff(values) != 1), length(values))

  hpdi <- lapply(seq(length(breaks) - 1), function(.x, values, breaks) {
    vals <- values[(breaks[.x] + 1):breaks[.x + 1]]
    vals <- unique(c(min(vals), max(vals)))
    paste0("[", paste(vals, collapse = ","), "]")
  }, values = values, breaks = breaks)

  paste(hpdi, collapse = ",")
}
