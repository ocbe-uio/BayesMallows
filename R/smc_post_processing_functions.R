#' @importFrom graphics mtext par

# Analysis

#' @title SMC Processing
#' @author Anja Stein
#' @param output input
#' @param colnames colnames
#' @importFrom gtools mixedorder
# AS: edited this function to include parameter `colnames`. This resolve issues in #118 with post processing functions not printing the names of items in rankings.
# The `default` is set to NULL so tat we do not cause plotting issues in `plot_rho_heatplot.
smc_processing<- function(output, colnames = NULL) {

  df <- data.frame(data = output)

  # if colnames are specified, then incorporate them
  if(is.null(colnames)){
    n_items <- ncol(df)
    cletters <- rep(c("Item"), times = n_items)
    cindexes <- (c(1:n_items))
    cnames <- c(paste(cletters, cindexes, sep = " "))
    colnames(df) <- cnames
  } else {
    colnames(df) <- colnames
  }

  new_df <- tidyr::gather(df, key = "item", value = "value")
  return(new_df)
}

heatmatrix <- function(output, burnin, rho){
  # get dimensions
  n_items <- dim(output)[2]
  sample <- (burnin+1):dim(output)[1]

  # get matrix of posterior probabilites by counting how many times an item is assigned specific rank
  heatmat <- matrix(nrow = n_items, ncol = n_items)
  colnames(output) <- NULL
  for (j in 1:n_items) {
    idx <- which(rho == j)
    heatmat[j, ] <- sapply(1:n_items, function(x) sum(output[sample, idx] == x))
  }

  # normalise to get probabilities
  heatmat <- heatmat / length(sample)
  return(heatmat)
}

# AS: adjusted font size values (`cex.lab` and 'cex`) of axis labels in heatmap function
create_heatplot <- function(mat, rho) {

  # sort items
  n_items <- length(rho)
  if (is.character(names(rho))) {
    items <- names(sort(rho))
  } else {
    items <- paste("0", order(rho), sep = "")
  }

  # create the heatplot
  par(mfrow = c(1, 1))
  fields::image.plot(mat,
    col = fields::tim.colors(64 * 10), axes = F, cex.lab = 1.0, zlim = range(0, 1),
    axis.args = list(at = seq(0, 1, 0.1), labels = seq(0, 1, 0.1), cex.axis = 1.2),
    xlab = "Consensus Ranking", ylab = "Rank"
  )

  # Add x and y axis labels and titles
  mtext(
    text = items, side = 1, line = 0.3,
    at = seq(0, 1, 1 / (n_items - 1)), las = 2, cex = 0.6
  )
  mtext(
    text = c(1:n_items), side = 2, line = 0.3,
    at = seq(0, 1, 1 / (n_items - 1)), las = 2, cex = 0.6
  )
  mtext(
    text = "Posterior probabilties for rho", side = 3, line = 0.3,
    at = , las = 1, cex = 1.2
  )
}

#' @title Plot rho heatplot
#' @param output input
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}
#' @param n_items Integer is the number of items in the consensus ranking
#' @param rho Numeric vector specifying the current consensus ranking
#' @return A heatplot
#' @author Anja Stein
#' @export
plot_rho_heatplot <- function(output, nmc, burnin, n_items, rho) {
  smc_plot <- smc_processing(output = output)

  # create the heatplot
  smc_rho_matrix <- matrix(smc_plot$value, ncol = n_items, nrow = nmc, byrow = FALSE)
  smc_heatmat_rho <- heatmatrix(output = smc_rho_matrix, burnin = burnin, rho = rho)
  smc_heatplot_full <- create_heatplot(mat = smc_heatmat_rho, rho = rho)
}

# same as compute_posterior_intervals, but removed the bayesmallows object and some other columns
compute_posterior_intervals_smc <- function(model_fit, burnin = model_fit$burnin,
                                               parameter = "alpha", level = 0.95,
                                               decimals = 3L) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)
  stopifnot(parameter %in% c("alpha", "rho", "cluster_probs", "cluster_assignment"))
  stopifnot(level > 0 && level < 1)


  if (burnin != 0) {
    df <- dplyr::filter(model_fit, .data$iteration > burnin) # removed model_fit[[parameter]]
  } else {
    df <- model_fit
  }

  if (parameter == "alpha" || parameter == "cluster_probs") {
    df <- dplyr::group_by(df, .data$cluster)
    df <- .compute_posterior_intervals(df, parameter, level, decimals)
  } else if (parameter == "rho") {
    decimals <- 0
    df <- dplyr::group_by(df, .data$cluster, .data$item)
    df <- .compute_posterior_intervals(df, parameter, level, decimals, discrete = TRUE)
  }

  df <- dplyr::ungroup(df)

  if (model_fit$n_clusters[1] == 1) df <- dplyr::select(df, -.data$cluster)

  return(df)
}


.compute_posterior_intervals <- function(df, parameter, level, decimals, discrete = FALSE) {
  dplyr::do(df, {
    format <- paste0("%.", decimals, "f")

    posterior_mean <- round(base::mean(.data$value), decimals)
    posterior_median <- round(stats::median(.data$value), decimals)

    if (discrete) {
      df <- dplyr::group_by(.data, .data$value)
      df <- dplyr::summarise(df, n = dplyr::n())
      df <- dplyr::arrange(df, dplyr::desc(.data$n))
      df <- dplyr::mutate(df,
        cumprob = cumsum(.data$n) / sum(.data$n),
        lagcumprob = dplyr::lag(.data$cumprob, default = 0)
      )

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
      if (is.matrix(hpdi)) {
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


#' Compute Consensus Ranking
#'
#' Compute the consensus ranking using either cumulative probability (CP) or maximum a posteriori (MAP) consensus
#' \insertCite{vitelli2018}{BayesMallows}. For mixture models, the
#' consensus is given for each mixture.
#'
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#'
#' @param type Character string specifying which consensus to compute. Either
#' \code{"CP"} or \code{"MAP"}. Defaults to \code{"CP"}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#' @author Anja Stein
#'
compute_consensus_smc <- function(model_fit, type, burnin) {
  if (type == "CP") {
    .compute_cp_consensus_smc(model_fit, burnin = burnin)
  } else if (type == "MAP") {
    .compute_map_consensus_smc(model_fit, burnin = burnin)
  }
}

.compute_cp_consensus_smc <- function(model_fit, burnin){
# FIXME: # 69 this function already exists on compute_consensus.R. Add S3 method?

  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)

  # Filter out the pre-burnin iterations

  if(burnin!=0){
    df <- dplyr::filter(model_fit, .data$iteration > burnin)
  }else {df <- model_fit}

  # Find the problem dimensions
  n_rows <- nrow(dplyr::distinct(df, .data$item, .data$cluster))

  # Check that there are rows.
  stopifnot(n_rows > 0)

  # Check that the number of rows are consistent with the information in
  # the model object
  stopifnot(model_fit$n_clusters * model_fit$n_items == n_rows)

  # Convert items and clustr to character, since factor levels are not needed in this case
  df <- dplyr::mutate_at(df, dplyr::vars(.data$item, .data$cluster), as.character)

  # Group by item, cluster, and value
  df <- dplyr::group_by(df, .data$item, .data$cluster, .data$value)

  # Find the count of each unique combination (value, item, cluster)
  df <- dplyr::count(df)

  # Arrange according to value, per item and cluster
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, .data$item, .data$cluster)
  df <- dplyr::arrange(df, .data$value, .by_group = TRUE)

  # Find the cumulative probability, by dividing by the total
  # count in (item, cluster) and the summing cumulatively
  df <- dplyr::mutate(df, cumprob = cumsum(.data$n/sum(.data$n)))

  # Find the CP consensus per cluster, using the find_cpc_smc function
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, .data$cluster)
  df <- dplyr::do(df, find_cpc_smc(.data))
  df <- dplyr::ungroup(df)

  # If there is only one cluster, we drop the cluster column
  if (model_fit$n_clusters[1] == 1) {
    df <- dplyr::select(df, -.data$cluster)
  }

  return(df)

}


# Internal function for finding CP consensus.
find_cpc_smc <- function(group_df){
# FIXME: # 69 this function already exists on compute_consensus.R. Add S3 method?
  # Declare the result dataframe before adding rows to it
  result <- dplyr::tibble(
    cluster = character(),
    ranking = numeric(),
    item = character(),
    cumprob = numeric()
  )
  n_items <- max(group_df$value)
  for(i in seq(from = 1, to = n_items, by = 1)){
    # Filter out the relevant rows
    tmp_df <- dplyr::filter(group_df, group_df$value == i)

    # Remove items in result
    tmp_df <- dplyr::anti_join(tmp_df, result, by = c("cluster", "item"))

    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    if (nrow(tmp_df) >= 1) {
      tmp_df <- dplyr::filter(tmp_df, .data$cumprob == max(.data$cumprob))
    }

    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)

    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(
      result,
      dplyr::select(
        tmp_df, .data$cluster, .data$ranking, .data$item, .data$cumprob
      )
    )

  }
  return(result)
}

 #AS: added one extra line of code to resolve of the issues in #118 with plotting too many rows in compute_rho_consensus
.compute_map_consensus_smc <- function(model_fit, burnin = model_fit$burnin){
# FIXME: # 69 this function already exists on compute_consensus.R. Add S3 method?

  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }

  if(burnin != 0){
    df <- dplyr::filter(model_fit, .data$iteration > burnin)
  } else {
    df <- model_fit
  }

  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))

  #-----------------------------------------------------------
  #AS: remove the column n_clusters, parameter
  df <- within(df, {n_clusters <- NULL; parameter <- NULL})
  #------------------------------------------------------------

  # Spread to get items along columns
  df <- tidyr::spread(df, key = .data$item, value = .data$value)

  # Group by everything except iteration, and count the unique combinations
  df <- dplyr::group_by_at(df, .vars = dplyr::vars(-.data$iteration))
  df <- dplyr::count(df)
  df <- dplyr::ungroup(df)
  # Keep only the maximum per cluster
  df <- dplyr::group_by(df, .data$cluster)
  df <- dplyr::mutate(df, n_max = max(.data$n))
  df <- dplyr::filter(df, .data$n == .data$n_max)
  df <- dplyr::ungroup(df)

  # Compute the probability
  df <- dplyr::mutate(df, probability = .data$n / n_samples)
  df <- dplyr::select(df, -.data$n_max, -.data$n)

  # Now collect one set of ranks per cluster
  df <- tidyr::gather(df, key = "item", value = "map_ranking",
                      -.data$cluster, -.data$probability)

  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, .data$cluster, .data$map_ranking)

  if (model_fit$n_clusters[1] == 1) {
    df <- dplyr::select(df, -.data$cluster)
  }

  return(df)

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
compute_posterior_intervals_rho <- function(output, nmc, burnin, colnames = NULL, verbose=FALSE) {
  #----------------------------------------------------------------
  # AS: added extra input parameter
  smc_plot <- smc_processing(output = output, colnames = colnames)
  #----------------------------------------------------------------
  smc_plot$n_clusters <- 1
  smc_plot$cluster <- "Cluster 1"

  rho_posterior_interval <- compute_posterior_intervals_smc(
    model_fit = smc_plot, burnin = burnin,
    parameter = "rho", level = 0.95, decimals = 2
  )

  #------------------------------------------------------------------------------------------
  #AS: reorder items to be in numerical order if no colnames are specified
  if (is.null(colnames)) {
    rho_posterior_interval  <- rho_posterior_interval[mixedorder(rho_posterior_interval$item), ]
  }
  #------------------------------------------------------------------------------------------

  if(verbose) print(rho_posterior_interval)
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
compute_rho_consensus <- function(output, nmc, burnin, C, type, colnames = NULL, verbose=FALSE) {

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

  # rho estimation using cumulative probability
  if (type == "CP") {
    results <- compute_consensus_smc(
      model_fit = smc_plot, type = "CP", burnin = burnin
    )
  } else {
    results <- compute_consensus_smc(
      model_fit = smc_plot, type = "MAP", burnin = burnin
    )
  }
  if (verbose) print(results)

  return(results)
}

#' @title Plot Alpha Posterior
#' @description posterior for alpha
#' @inheritParams compute_posterior_intervals_rho
#' @export
#' @author Anja Stein
#'
# AS: if you remove the verbose input variable, then the function will be consistent
# with the other plot functions(they all print when verbose=FALSE, but this function doesn't.)
#`plot_rho_heatplot` doesn't require the variable `verbose`,
# so I'm not sure if this function does to plot the density of alpha
plot_alpha_posterior <- function(output, nmc, burnin) {
  alpha_samples_table <- data.frame(iteration = 1:nmc, value = output)

  plot_posterior_alpha <- ggplot2::ggplot(alpha_samples_table, ggplot2::aes(x = alpha_samples_table$value)) +
    ggplot2::geom_density() +
    ggplot2::xlab(expression(alpha)) +
    ggplot2::ylab("Posterior density") +
    ggplot2::ggtitle(label = "Implemented SMC scheme") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

  print(plot_posterior_alpha)
}

#' @title Compute Posterior Intervals Alpha
#' @description posterior confidence intervals
#' @inheritParams compute_posterior_intervals_rho
#' @export
#' @author Anja Stein
#'
compute_posterior_intervals_alpha <- function(output, nmc, burnin, verbose=FALSE) {
  alpha_samples_table <- data.frame(iteration = 1:nmc, value = output)
  alpha_samples_table$n_clusters <- 1
  alpha_samples_table$cluster <- "Cluster 1"


  alpha_mixture_posterior_interval <- compute_posterior_intervals_smc(alpha_samples_table,
    burnin = burnin,
    parameter = "alpha", level = 0.95, decimals = 2
  )
  if (verbose) print(alpha_mixture_posterior_interval)
  return(alpha_mixture_posterior_interval)
}
                      
                      
                      
 # plot the posterior for rho for each item
#' @param output input
#' @param nmc Number of Monte Carlo samples
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}
#' @param C Number of cluster
#' @colnames A vector of item names. If NULL, we generate generic names for the items in the ranking.
#' @param items Either a vector of item names, or a
#'   vector of indices. If NULL, five items are selected randomly.
plot_rho_posterior <- function(output, nmc, burnin, C, colnames = NULL, items = NULL){

  n_items = dim(output)[2]
  
  if(is.null(items) && n_items > 5){
    message("Items not provided by user or more than 5 items in a ranking. Picking 5 at random.")
    items <- sample(1:n_items, 5, replace = F)
    items = sort(items)
    
  } else if (is.null(items) && n_items <= 5) {
    items <- c(1:n_items)
    items = sort(items)
  }
  
  # do smc processing here
  smc_plot = smc_processing(output = output, colnames = colnames)
  
  if(!is.character(items)){
    items <- unique(smc_plot$item)[items]
  }
  
  iteration = rep(c(1:nmc), times = n_items)
  df = cbind(iteration, smc_plot)
  
  if(C==1){
    df = cbind(cluster = "Cluster 1", df)
  }
   
  df <- dplyr::filter(df, .data$iteration > burnin, .data$item %in% items)
  
  # Compute the density, rather than the count, since the latter
  # depends on the number of Monte Carlo samples
  df <- dplyr::group_by(df, .data$cluster, .data$item, .data$value)
  df <- dplyr::summarise(df, n = dplyr::n())
  df <- dplyr::mutate(df, pct = .data$n / sum(.data$n))
  
  df$item <- factor(df$item, levels = c(items))
  
  # Taken from misc.R function in BayesMallows
  scalefun <- function(x) sprintf("%d", as.integer(x))
  
  # Finally create the plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$pct)) +
    ggplot2::geom_col() +
    ggplot2::scale_x_continuous(labels = scalefun) +
    ggplot2::xlab("rank") +
    ggplot2::ylab("Posterior probability")
  
  if(C == 1){
    p <- p + ggplot2::facet_wrap(~ .data$item)
  } else {
    p <- p + ggplot2::facet_wrap(~ .data$cluster + .data$item)
  }
  
  return(p)                     
                     
