#
#' Compute CP Consensus
#'
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#' @param burnin The number of iterations to discard as burn-in.
#'
#' @details Computes the cumulative probability (CP) consensus rankings of the latent ranks.
#' If the number of clusters is larger than 1, it returns the CP consensus per cluster.
#'
#' @export
#'
compute_cp_consensus <- function(model_fit, burnin){
  stopifnot(class(model_fit) == "BayesMallows")

  # Filter out the pre-burnin iterations
  df <- dplyr::filter(model_fit$rho, .data$iteration > burnin)

  # Find the problem dimensions
  n_rows <- nrow(dplyr::distinct(df, .data$item, .data$cluster))

  # Check that there are rows.
  stopifnot(n_rows > 0)

  # Check that the number of rows are consistent with the information in
  # the model object
  stopifnot(model_fit$n_clusters * model_fit$n_items == n_rows)

  # Group by item, cluster, and value
  df <- dplyr::group_by(df, .data$item, .data$cluster, .data$value)

  # Find the count of each unique combination (value, item, cluster)
  df <- dplyr::count(df)

  # Divide by the number of counts in total, per item and cluster
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, .data$item, .data$cluster)
  df <- dplyr::mutate(df, prob = .data$n/sum(.data$n))

  # Arrange according to value, per item and cluster
  df <- dplyr::arrange(df, .data$value, .by_group = TRUE)

  # Compute the cumulative probability per item and cluster
  df <- dplyr::mutate(df,  cumprob = cumsum(.data$prob))
  df <- dplyr::ungroup(df)

  # We now have reduced df to (n_items * n_clusters) rows, which
  # is typically a small number.
  # Declare the result dataframe before adding rows to it
  result <- dplyr::tibble()

  for(i in seq(from = 1, to = model_fit$n_items, by = 1)){
    # Filter out the relevant rows.
    tmp_df <- dplyr::filter(df, .data$value == i)

    # Keep the max per cluster only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    tmp_df <- dplyr::group_by(tmp_df, .data$cluster)
    tmp_df <- dplyr::filter(tmp_df, .data$cumprob == max(.data$cumprob))

    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)

    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(result,
                               dplyr::select(tmp_df, .data$cluster, .data$ranking, .data$item, .data$cumprob))

    # Remove the select rows from df
    # Here we remove everything ranked i
    df <- dplyr::filter(df, .data$value != i)
    # Here we remove the (item, cluster) that are in tmp_df
    df <- dplyr::anti_join(df, tmp_df, by = c("cluster", "item"))
  }

  # If there is only one cluster, we drop the cluster column
  if(model_fit$n_clusters == 1){
    result <- dplyr::select(result, -.data$cluster)
  }

  return(result)

}
