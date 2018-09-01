#' Maximum a Posterior Consensus Ranking
#'
#' Compute the maximum a posterior consensus ranking. If a mixture model
#' is supplied, the maximum a posterior consensus ranking for each mixture
#' is returned.
#'
#' @param model_fit An object of type \code{BayesMallows}, returned from
#'   \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. See \code{\link{assess_convergence}}.
#'
#' @export
#'
#' @example /inst/examples/compute_map_consensus_example.R
#'
compute_map_consensus <- function(model_fit, burnin){
  df <- dplyr::filter(model_fit$rho, .data$iteration > burnin)

  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))

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

  if(model_fit$n_clusters == 1){
    df <- dplyr::select(df, -.data$cluster)
  }

  return(df)

}
