#' Compute CP Consensus
#'
#' Compute the cumulative probability (CP) consensus ranking
#' \insertCite{vitelli2018}{BayesMallows} of the items. For mixture models, the
#' CP consensus is given for each ranking.
#'
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#'
#' @param burnin A numeric value specifying the number of iterations to discard
#'   as burn-in. See \code{\link{assess_convergence}}.
#'
#' @references \insertAllCited{}
#'
#'
#' @export
#'
#' @examples
#' # Analysis of complete rankings
#' # The example datasets potato_visual and potato_weighing contain complete
#' # rankings of 20 items, by 12 assessors. We first analyse these using the Mallows
#' # model:
#' model_fit <- compute_mallows(potato_visual)
#' # We study the trace plot of the parameters
#' # alpha is the default
#' assess_convergence(model_fit)
#' # When studying convergence of rho, we can also specify which items to plot
#' assess_convergence(model_fit, type = "rho", items = 1:5)
#' # Based on these plots, we conclude that the Markov chain has converged well
#' # before 1,000 iterations. We hence set burnin = 1000.
#' # Next, we use the generic plot function to study the posterior distributions
#' # of alpha and rho
#' plot(model_fit, burnin = 1000)
#' plot(model_fit, burnin = 1000, type = "rho", items = 1:20)
#' # We can also compute the CP consensus posterior ranking
#' compute_cp_consensus(model_fit, burnin = 1000)
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

  # Convert items to character, since factor levels are not needed in this case
  df <- dplyr::mutate(df, item = as.character(.data$item))

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

  # Find the CP consensus per cluster, using the find_cpc function
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, .data$cluster)
  df <- dplyr::do(df, find_cpc(.data))
  df <- dplyr::ungroup(df)

  # If there is only one cluster, we drop the cluster column
  if(model_fit$n_clusters == 1){
    df <- dplyr::select(df, -.data$cluster)
  }

  return(df)

}


# Internal function for finding CP consensus.
find_cpc <- function(group_df){
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
    tmp_df <- dplyr::filter(group_df, .data$value == i)

    # Remove items in result
    tmp_df <- dplyr::anti_join(tmp_df, result, by = c("cluster", "item"))

    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    tmp_df <- dplyr::filter(tmp_df, .data$cumprob == max(.data$cumprob))

    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)

    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(result,
                               dplyr::select(tmp_df, .data$cluster, .data$ranking, .data$item, .data$cumprob))

  }
  return(result)
}
