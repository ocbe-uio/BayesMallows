#' Compute Consensus Ranking
#'
#' Compute the consensus ranking using either cumulative probability (CP) or maximum a posteriori (MAP) consensus
#' \insertCite{vitelli2018}{BayesMallows}. For mixture models, the
#' consensus is given for each mixture. Consensus of augmented ranks can also be computed
#' for each assessor, by setting \code{parameter = "Rtilde"}.
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
#' @param parameter Character string defining the parameter for which to compute the
#' consensus. Defaults to \code{"rho"}. Available options are \code{"rho"} and \code{"Rtilde"},
#' with the latter giving consensus rankings for augmented ranks.
#'
#' @param assessors When \code{parameter = "rho"}, this integer vector is used to
#' define the assessors for which to compute the augmented ranking. Defaults to
#' \code{1L}, which yields augmented rankings for assessor 1.
#'
#' @references \insertAllCited{}
#'
#'
#' @export
#'
#' @example /inst/examples/compute_consensus_example.R
#'
compute_consensus <- function(model_fit, type = "CP", burnin = model_fit$burnin,
                              parameter = "rho", assessors = 1L){

  if(is.null(burnin)){
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  stopifnot(class(model_fit) == "BayesMallows")

  if(parameter == "Rtilde" && !inherits(model_fit$augmented_data, "data.frame")){
    stop(paste("For augmented ranks, please refit",
               "model with option 'save_aug = TRUE'."))
  }

  if(parameter == "rho"){
    # Filter out the pre-burnin iterations
    df <- dplyr::filter(model_fit$rho, .data$iteration > burnin)

    # Find the problem dimensions
    n_rows <- nrow(dplyr::distinct(df, .data$item, .data$cluster))

    # Check that there are rows.
    stopifnot(n_rows > 0)

    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(model_fit$n_clusters * model_fit$n_items == n_rows)

    df <- if(type == "CP"){
      .compute_cp_consensus(df)
    } else if(type == "MAP"){
      .compute_map_consensus(df)
    }


  } else if(parameter == "Rtilde"){
    # Filter out the pre-burnin iterations and get the right assessors
    df <- dplyr::filter(model_fit$augmented_data, .data$iteration > burnin, .data$assessor %in% assessors)

    # Find the problem dimensions
    n_rows <- nrow(dplyr::distinct(df, .data$assessor, .data$item))

    # Check that there are rows.
    stopifnot(n_rows > 0)

    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(length(assessors) * model_fit$n_items == n_rows)

    # Treat assessors as clusters
    df <- dplyr::rename(df, cluster = "assessor")

    df <- if(type == "CP"){
      .compute_cp_consensus(df)
    } else if(type == "MAP"){
      .compute_map_consensus(df)
    }

    if("cluster" %in% names(df)){
      df <- dplyr::rename(df, assessor = "cluster")
    }

  }

  return(df)


}

.compute_cp_consensus <- function(df){

  # Convert items and cluster to character, since factor levels are not needed in this case
  df <- dplyr::mutate_at(df, dplyr::vars(.data$item, .data$cluster),
                         as.character)

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
  if(length(unique(df$cluster)) == 1){
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
    tmp_df <- dplyr::filter(tmp_df, .data$cumprob == max(dplyr::coalesce(.data$cumprob, 0)))

    # Add the ranking
    tmp_df <- dplyr::mutate(tmp_df, ranking = i)

    # Select the columns we want to keep, and put them in result
    result <- dplyr::bind_rows(result,
                               dplyr::select(tmp_df, .data$cluster, .data$ranking, .data$item, .data$cumprob))

  }
  return(result)
}

.compute_map_consensus <- function(df){

  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))

  # Spread to get items along columns
  df <- tidyr::pivot_wider(df, names_from = "item", values_from = "value")

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
  df <- tidyr::pivot_longer(df, cols = -c("cluster", "probability"), names_to = "item", values_to = "map_ranking")

  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, .data$cluster, .data$map_ranking)

  if(length(unique(df$cluster)) == 1){
    df <- dplyr::select(df, -.data$cluster)
  }



  return(df)

}
