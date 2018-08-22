


#' Compute within-cluster distance to latent ranking
#'
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#' @param rankings A matrix with complete rankings, the same as was
#' given to \code{\link{compute_mallows}}.
#' @param burnin The number of iterations to discard.
#' @param probs Probability of the quantiles to returned. This argument is passed on to
#' \code{stats::quantile}.
#' @param thinning Optional argument which can be used to take only every
#' \code{thinning}th iteration after burn-in. If the function takes very long to run,
#' one might try to set this number to something larger than 1. Be aware that this always leads to
#' less accurate posterior distributions.
#'
#' @details At the moment, this function computes the footrule distance,
#' and no other distance measures are supported.
#'
#' @return A vector of quantiles.
#' @export
#'
compute_within_cluster_distance <- function(model_fit, rankings, burnin,
                                            probs = seq(0, 1, 0.25),
                                            thinning = 1){

  df <- .compute_within_cluster_distance(model_fit, rankings, burnin, thinning)
  return(stats::quantile(df$within_cluster_distance, probs = probs))

}


# Internal function for computing with cluster distance
# Defined this way because it is used by more than one function
.compute_within_cluster_distance <- function(model_fit, rankings, burnin, thinning){
  stopifnot(burnin < model_fit$nmc)

  # Extract the values and filter out pre-burnin iterations
  cluster_indicator <- model_fit$cluster_indicator
  cluster_indicator <- dplyr::filter(cluster_indicator, .data$iteration > burnin, .data$iteration %% thinning == 0)
  cluster_indicator <- dplyr::rename(cluster_indicator, cluster = .data$value)

  rho_c <- model_fit$rho
  rho_c <- dplyr::filter(rho_c, .data$iteration > burnin, .data$iteration %% thinning == 0)

  # Convert the rankings to "long" format
  rank_df <- dplyr::as_tibble(rankings)
  rank_df <- dplyr::mutate(rank_df, assessor = dplyr::row_number())
  rank_df <- tidyr::gather(rank_df, key = "item", value = "value", -.data$assessor)
  rank_df <- dplyr::rename(rank_df, ranking = .data$value)

  # Join the dataframes
  df <- dplyr::inner_join(cluster_indicator, rank_df, by = "assessor")
  df <- dplyr::inner_join(df, rho_c, by = c("cluster", "item", "iteration"))

  # Find the footrule distance per iteration
  df <- dplyr::group_by(df, .data$iteration)
  df <- dplyr::summarise(df, within_cluster_distance = sum(abs(.data$ranking - .data$value)))

  # Add the number of clusters
  df <- dplyr::mutate(df, n_clusters = model_fit$n_clusters)

  return(df)
}

