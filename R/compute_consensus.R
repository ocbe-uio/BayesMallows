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
#'
#' @references \insertAllCited{}
#'
#'
#' @export
#'
#' @example /inst/examples/compute_consensus_example.R
#'
compute_consensus <- function(model_fit, type = "CP", burnin = model_fit$burnin){

  if(type == "CP"){
    compute_cp_consensus(model_fit, burnin = burnin)
  } else if(type == "MAP"){
    compute_map_consensus(model_fit, burnin = burnin)
  }


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
