#' @importFrom graphics mtext par

#' @title SMC Processing
#' @author Anja Stein
#' @param output input
#' @param colnames colnames
# AS: edited this function to include parameter `colnames`. This resolve issues in #118 with post processing functions not printing the names of items in rankings.
# The `default` is set to NULL so tat we do not cause plotting issues in `plot_rho_heatplot.
smc_processing <- function(output, colnames = NULL) {
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
  new_df$id <- NULL # drop the "id" column
  class(new_df) <- c("SMCMallows", "data.frame")
  return(new_df)
}


#' @title Compute Posterior Intervals Rho
#' @description posterior confidence intervals for rho
#' @inheritParams smc_processing
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#' @param verbose if \code{TRUE}, prints the final output even if the function
#' is assigned to an object. Defaults to \code{FALSE}.
#' @export
#' @author Anja Stein
#'
# AS: added an extra inout variable `colnames`. This is called in the function `smc_processing`.
compute_posterior_intervals_rho <- function(output, nmc, burnin, colnames = NULL, verbose = FALSE) {
  #----------------------------------------------------------------
  # AS: added extra input parameter
  smc_plot <- smc_processing(output = output, colnames = colnames)
  #----------------------------------------------------------------
  smc_plot$n_clusters <- 1
  smc_plot$cluster <- "Cluster 1"

  rho_posterior_interval <- compute_posterior_intervals(
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
#' @description MAP AND CP consensus ranking estimates
#' @inheritParams compute_posterior_intervals_rho
#' @param C C
#' @param type type
#' @export
#' @author Anja Stein
#'
# AS: added an extra inout variable `colnames`. This is called in the function `smc_processing`.
compute_rho_consensus <- function(output, nmc, burnin, C, type, colnames = NULL, verbose = FALSE) {
  n_items <- dim(output)[2]

  #----------------------------------------------------------------
  # AS: added extra input parameter
  smc_plot <- smc_processing(output = output, colnames = colnames)
  #----------------------------------------------------------------

  iteration <- array(rep((1:nmc), n_items))
  smc_plot <- data.frame(data = cbind(iteration, smc_plot))
  colnames(smc_plot) <- c("iteration", "item", "value")

  smc_plot$n_clusters <- C
  smc_plot$parameter <- "rho"
  smc_plot$cluster <- "cluster 1"

  class(smc_plot) <- c("consensus_SMCMallows", "data.frame")

  # rho estimation using cumulative probability
  if (type == "CP") {
    results <- compute_consensus(
      model_fit = smc_plot, type = "CP", burnin = burnin
    )
  } else {
    results <- compute_consensus(
      model_fit = smc_plot, type = "MAP", burnin = burnin
    )
  }
  if (verbose) print(results)

  return(results)
}

#' @title Compute Posterior Intervals Alpha
#' @description posterior confidence intervals
#' @inheritParams compute_posterior_intervals_rho
#' @export
#' @author Anja Stein
#'
compute_posterior_intervals_alpha <- function(output, nmc, burnin, verbose = FALSE) {
  alpha_samples_table <- data.frame(iteration = 1:nmc, value = output)
  alpha_samples_table$n_clusters <- 1
  alpha_samples_table$cluster <- "Cluster 1"
  class(alpha_samples_table) <- c("SMCMallows", "data.frame")

  alpha_mixture_posterior_interval <- compute_posterior_intervals(alpha_samples_table,
    burnin = burnin,
    parameter = "alpha", level = 0.95, decimals = 2
  )
  if (verbose) print(alpha_mixture_posterior_interval)
  return(alpha_mixture_posterior_interval)
}
