#' @title SMC Processing
#' @param output a subset of an SMCMallows object (though technically any matrix will do)
#' @param colnames colnames
#' @return A processed file of the SMCMallows class
#' @seealso \code{\link{smc_mallows_new_item_rank}} and
#' \code{\link{smc_mallows_new_users}}, which are functions generating objects
#' of SMCMallows class.
#' @noRd
#' @importFrom methods is
smc_processing <- function(output, colnames = NULL) {
  # Recasting of input for proper handling below
  df <- data.frame(data = output)

  # if colnames are specified, then incorporate them
  if (is.null(colnames)) {
    n_items <- ncol(df)
    cletters <- rep("Item", times = n_items)
    cindexes <- seq_len(n_items)
    cnames <- c(paste(cletters, cindexes, sep = " "))
    colnames(df) <- cnames
  } else {
    colnames(df) <- colnames
  }
  new_df <- stats::reshape(
    df,
    direction = "long",
    varying = names(df),
    new.row.names = seq_len(prod(dim(df))),
    v.names = "value",
    timevar = "item",
    times = names(df)
  )
  new_df$id <- NULL # the "id" should not be part of the SMCMallows object
  class(new_df) <- c("SMCMallows", "data.frame")
  return(new_df)
}

#' @title Compute Posterior Intervals Rho
#' @description posterior confidence intervals for rho
#' @param output a subset of an SMCMallows object (though technically any matrix
#'   will do)
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in.
#' @param colnames Column names.
#' @param verbose if \code{TRUE}, prints the final output even if the function
#' is assigned to an object. Defaults to \code{FALSE}.
#' @export
#' @inherit compute_rho_consensus examples
#' @family deprecated
compute_posterior_intervals_rho <- function(output, nmc, burnin, colnames = NULL, verbose = FALSE) {
  .Deprecated(
    new = "compute_posterior_intervals",
    msg = paste(
      "compute_posterior_intervals_rho() is deprecated. Please",
      "use compute_posterior_intervals() with argument parameter = 'rho'."
    )
  )
  # Validation
  stopifnot(is(output, "matrix"))

  #----------------------------------------------------------------
  # AS: added extra input parameter
  smc_plot <- smc_processing(output = output, colnames = colnames)
  #----------------------------------------------------------------
  smc_plot$n_clusters <- 1
  smc_plot$cluster <- "Cluster 1"

  rho_posterior_interval <- compute_posterior_intervals_SMCMallows_deprecated(
    model_fit = smc_plot, burnin = burnin,
    parameter = "rho", level = 0.95, decimals = 2
  )

  #------------------------------------------------------------------------------------------
  # AS: reorder items to be in numerical order if no colnames are specified
  if (is.null(colnames)) {
    item_numbers <- as.numeric(gsub("\\D", "", rho_posterior_interval$item))
    mixed_order <- match(sort(item_numbers), item_numbers)
    rho_posterior_interval <- rho_posterior_interval[mixed_order, ]
  }
  #------------------------------------------------------------------------------------------

  if (verbose) print(rho_posterior_interval)
  return(rho_posterior_interval)
}

#' @title Compute rho consensus
#'
#' @description This function has been deprecated. Please use
#'   \code{\link{compute_consensus.SMCMallows}}.
#'
#'
#' @description MAP and CP consensus ranking estimates
#'
#'   This function is deprecated. Please use
#'   \code{\link{compute_consensus.SMCMallows}} instead.
#'
#' @param C C Number of clusters.
#' @param output Matrix
#' @param nmc Number of Monte Carlo samples
#' @param burnin Burnin
#' @param colnames Column names
#' @param verbose Logical
#' @param type type
#'
#' @export
#' @family deprecated
compute_rho_consensus <- function(output, nmc, burnin, C, type = "CP", colnames = NULL, verbose = FALSE) {
  .Deprecated("compute_consensus",
    msg = paste(
      "compute_rho_consensus has been deprecated. Please use",
      "compute_consensus instead."
    )
  )
}

#' @title Compute Posterior Intervals Alpha
#' @description posterior confidence intervals
#' @inheritParams compute_posterior_intervals_rho
#' @export
#' @inherit compute_rho_consensus examples
#' @family deprecated
compute_posterior_intervals_alpha <- function(output, nmc, burnin, verbose = FALSE) {
  .Deprecated(
    new = "compute_posterior_intervals",
    msg = paste(
      "compute_posterior_intervals_alpha() is deprecated. Please",
      "use compute_posterior_intervals() with argument parameter = 'alpha'."
    )
  )
  # Validation
  stopifnot(is(output, "numeric"))

  alpha_samples_table <- data.frame(iteration = 1:nmc, value = output)
  alpha_samples_table$n_clusters <- 1
  alpha_samples_table$cluster <- "Cluster 1"
  class(alpha_samples_table) <- c("SMCMallows", "data.frame")

  alpha_mixture_posterior_interval <- compute_posterior_intervals_SMCMallows_deprecated(alpha_samples_table,
    burnin = burnin,
    parameter = "alpha", level = 0.95, decimals = 2
  )
  if (verbose) print(alpha_mixture_posterior_interval)
  return(alpha_mixture_posterior_interval)
}


compute_posterior_intervals_SMCMallows_deprecated <- function(model_fit, burnin = model_fit$burnin, parameter = "alpha", level = 0.95,
                                                              decimals = 3L, ...) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)
  stopifnot(parameter %in% c("alpha", "rho", "cluster_probs", "cluster_assignment"))
  stopifnot(level > 0 && level < 1)


  if (burnin != 0) {
    df <- model_fit[model_fit$iteration > burnin, , drop = FALSE]
  } else {
    df <- model_fit
  }

  if (parameter == "alpha" || parameter == "cluster_probs") {
    df <- .compute_posterior_intervals(split(df, f = df$cluster), parameter, level, decimals)
  } else if (parameter == "rho") {
    decimals <- 0
    df <- .compute_posterior_intervals(split(df, f = interaction(df$cluster, df$item)),
      parameter, level, decimals,
      discrete = TRUE
    )
  }


  if (model_fit$n_clusters[1] == 1) df$cluster <- NULL

  return(df)
}
