#' @title Compute Consensus Ranking
#' @description Compute the consensus ranking using either cumulative
#' probability (CP) or maximum a posteriori (MAP) consensus
#' \insertCite{vitelli2018}{BayesMallows}. For mixture models, the
#' consensus is given for each mixture. Consensus of augmented ranks
#' can also be computed
#' for each assessor, by setting \code{parameter = "Rtilde"}.
#' @param model_fit An object returned from \code{\link{compute_mallows}}.
#' @param ... other arguments passed to methods.
#' @references \insertAllCited{}
#' @export
#' @example /inst/examples/compute_consensus_example.R
compute_consensus <- function(model_fit, ...) {
  UseMethod("compute_consensus")
}

.compute_cp_consensus <- function(df, ...) {
  UseMethod(".compute_cp_consensus")
}

.compute_map_consensus <- function(df, ...) {
  UseMethod(".compute_map_consensus")
}

find_cpc <- function(df, ...) {
  UseMethod("find_cpc")
}

#' @title Compute Consensus Ranking
#' @inheritParams compute_consensus
#' @param type Character string specifying which consensus to compute. Either
#' \code{"CP"} or \code{"MAP"}. Defaults to \code{"CP"}.
#' @param burnin A numeric value specifying the number of iterations
#' to discard as burn-in. Defaults to \code{model_fit$burnin}, and must be
#' provided if \code{model_fit$burnin} does not exist. See \code{\link{assess_convergence}}.
#' @param parameter Character string defining the parameter for which to compute the
#' consensus. Defaults to \code{"rho"}. Available options are \code{"rho"} and \code{"Rtilde"},
#' with the latter giving consensus rankings for augmented ranks.
#' @param assessors When \code{parameter = "rho"}, this integer vector is used to
#' define the assessors for which to compute the augmented ranking. Defaults to
#' \code{1L}, which yields augmented rankings for assessor 1.
#' @export
compute_consensus.BayesMallows <- function(
  model_fit, type = "CP", burnin = model_fit$burnin, parameter = "rho",
  assessors = 1L, ...
) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }
  stopifnot(burnin < model_fit$nmc)

  stopifnot(inherits(model_fit, "BayesMallows"))

  if (parameter == "Rtilde" && !inherits(model_fit$augmented_data, "data.frame")) {
    stop("For augmented ranks, please refit model with option 'save_aug = TRUE'.")
  }

  if (parameter == "rho") {
    # Filter out the pre-burnin iterations
    df <- model_fit$rho[model_fit$rho$iteration > burnin, , drop = FALSE]

    # Find the problem dimensions
    n_rows <- length(unique(paste(df$item, df$cluster)))

    # Check that there are rows.
    stopifnot(n_rows > 0)

    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(model_fit$n_clusters * model_fit$n_items == n_rows)

    class(df) <- c("consensus_BayesMallows", "tbl_df", "tbl", "data.frame")

    df <- if (type == "CP") {
      .compute_cp_consensus(df)
    } else if (type == "MAP") {
      .compute_map_consensus(df)
    }


  } else if (parameter == "Rtilde") {
    # Filter out the pre-burnin iterations and get the right assessors
    df <- model_fit$augmented_data[model_fit$augmented_data$iteration > burnin &
                                     model_fit$augmented_data$assessor %in% assessors, , drop = FALSE]

    # Find the problem dimensions
    n_rows <- length(unique(paste(df$assessor, df$item)))

    # Check that there are rows.
    stopifnot(n_rows > 0)

    # Check that the number of rows are consistent with the information in
    # the model object
    stopifnot(length(assessors) * model_fit$n_items == n_rows)

    # Treat assessors as clusters
    names(df)[names(df) == "assessor"] <- "cluster"
    class(df) <- c("consensus_BayesMallows", "tbl_df", "tbl", "data.frame")

    df <- if (type == "CP") {
      .compute_cp_consensus(df)
    } else if (type == "MAP") {
      .compute_map_consensus(df)
    }

    if ("cluster" %in% names(df)) {
      names(df)[names(df) == "cluster"] <- "assessor"
    }

  }

  return(df)

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
compute_consensus.consensus_SMCMallows <- function(model_fit, type, burnin) {
  if (type == "CP") {
    .compute_cp_consensus(model_fit, burnin = burnin)
  } else if (type == "MAP") {
    .compute_map_consensus(model_fit, burnin = burnin)
  }
}

.compute_cp_consensus.consensus_BayesMallows <- function(df) {
  # Count per item, cluster, and value
  df <- aggregate(
    list(n = df$iteration),
    by = list(item = as.character(df$item),
              cluster = as.character(df$cluster), value = df$value),
    FUN = length)

  # Arrange according to value, per item and cluster
  df <- do.call(rbind, lapply(split(df, f = ~ item + cluster), function(x){
    x <- x[order(x$value), ]
    x$cumprob <- cumsum(x$n) / sum(x$n)
    x
  }))

  # Find the CP consensus per cluster, using the find_cpc function
  df <- dplyr::group_by(df, .data$cluster)
  class(df) <- c("consensus_BayesMallows", "grouped_df", "tbl_df", "tbl", "data.frame")
  df <- find_cpc(df)
  df <- dplyr::ungroup(df)

  df <- df[order(df$cluster, df$ranking), ]

  # If there is only one cluster, we drop the cluster column
  if (length(unique(df$cluster)) == 1) {
    df$cluster <- NULL
  }

  return(df)
}

.compute_cp_consensus.consensus_SMCMallows <- function(model_fit, burnin) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }

  stopifnot(burnin < model_fit$nmc)

  # Filter out the pre-burnin iterations

  if (burnin != 0) {
    df <- model_fit[model_fit$iteration > burnin, , drop = FALSE]
  } else {
    df <- model_fit
  }

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
  df <- dplyr::mutate(df, cumprob = cumsum(.data$n / sum(.data$n)))

  # Find the CP consensus per cluster, using the find_cpc function
  df <- dplyr::ungroup(df)
  df <- dplyr::group_by(df, .data$cluster)
  class(df) <- c("consensus_SMCMallows", "grouped_df", "tbl_df", "tbl", "data.frame")
  df <- find_cpc(df)
  df <- dplyr::ungroup(df)

  # If there is only one cluster, we drop the cluster column
  if (model_fit$n_clusters[1] == 1) {
    df$cluster <- NULL
  }

  return(df)

}

# Internal function for finding CP consensus.
find_cpc.consensus_BayesMallows <- function(group_df, group_var = "cluster") {
  # Declare the result dataframe before adding rows to it
  result <- data.frame(
    cluster = character(),
    ranking = numeric(),
    item = character(),
    cumprob = numeric()
  )
  n_items <- max(group_df$value)
  group_df$cumprob[is.na(group_df$cumprob)] <- 0

  for (i in seq(from = 1, to = n_items, by = 1)) {
    # Filter out the relevant rows
    tmp_df <- group_df[group_df$value == i, , drop = FALSE]

    # Remove items in result
    tmp_df <- tmp_df[!interaction(tmp_df[c("cluster", "item")]) %in%
                       interaction(result[c("cluster", "item")]), ]

    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    tmp_df <- do.call(rbind,
                      lapply(split(tmp_df, f = tmp_df[group_var]), function(x) {
                        x[x$cumprob == max(x$cumprob), ]} ))

    # Add the ranking
    tmp_df$ranking <- i

    # Select the columns we want to keep, and put them in result
    result <- rbind(
      result,
      tmp_df[, c("cluster", "ranking", "item", "cumprob"), drop = FALSE])

  }
  return(result)
}

# Internal function for finding CP consensus.
find_cpc.consensus_SMCMallows <- function(group_df) {
  # Declare the result dataframe before adding rows to it
  result <- dplyr::tibble(
    cluster = character(),
    ranking = numeric(),
    item = character(),
    cumprob = numeric()
  )
  n_items <- max(group_df$value)
  for (i in seq(from = 1, to = n_items, by = 1)) {
    # Filter out the relevant rows
    tmp_df <- group_df[group_df$value == i, , drop = FALSE]

    # Remove items in result
    tmp_df <- dplyr::anti_join(tmp_df, result, by = c("cluster", "item"))

    # Keep the max only. This filtering must be done after the first filter,
    # since we take the maximum among the filtered values
    if (nrow(tmp_df) >= 1) {
      tmp_df <- tmp_df[tmp_df$cumprob == max(tmp_df$cumprob), , drop = FALSE]
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

.compute_map_consensus.consensus_BayesMallows <- function(df) {

  # Store the total number of iterations after burnin
  n_samples <- length(unique(df$iteration))

  # Reshape to get items along columns
  df <- stats::reshape(as.data.frame(df), direction = "wide",
                       idvar = c("cluster", "iteration"),
                       timevar = "item")
  names(df) <- gsub("^value\\.", "", names(df))

  # Group by everything except iteration, and count the unique combinations
  df <- dplyr::group_by_at(df, .vars = dplyr::vars(-.data$iteration))
  df <- dplyr::count(df)
  df <- dplyr::ungroup(df)
  # Keep only the maximum per cluster
  df <- dplyr::group_by(df, .data$cluster)
  df <- dplyr::mutate(df, n_max = max(.data$n))
  df <- df[df$n == df$n_max, , drop = FALSE]
  df <- dplyr::ungroup(df)

  # Compute the probability
  df$probability <- df$n / n_samples
  df$n_max <- NULL
  df$n <- NULL

  # Now collect one set of ranks per cluster
  df$id <- seq_len(nrow(df))
  df <- stats::reshape(as.data.frame(df), direction = "long",
          varying = setdiff(names(df), c("cluster", "probability", "id")),
          v.names = "map_ranking",
          timevar = "item",
          idvar = c("cluster", "probability", "id"),
          times = setdiff(names(df), c("cluster", "probability", "id")))
  rownames(df) <- NULL
  df$id <- NULL

  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, .data$cluster, .data$map_ranking)

  if (length(unique(df$cluster)) == 1) {
    df <- dplyr::select(df, -.data$cluster)
  }



  return(df)

}

 #AS: added one extra line of code to resolve of the issues in #118 with plotting too many rows in compute_rho_consensus
.compute_map_consensus.consensus_SMCMallows <- function(model_fit, burnin = model_fit$burnin) {
  if (is.null(burnin)) {
    stop("Please specify the burnin.")
  }

  if (burnin != 0) {
    df <- model_fit[model_fit$iteration > burnin, , drop = FALSE]
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
  df <- stats::reshape(
    data = as.data.frame(df),
    direction = "wide",
    idvar = c("iteration", "cluster"),
    timevar = "item",
    varying = list(unique(df$item))
  )
  attr(df, "reshapeWide") <- NULL # maintain identity to spread() output

  # Group by everything except iteration, and count the unique combinations
  df <- dplyr::group_by_at(df, .vars = dplyr::vars(-.data$iteration))
  df <- dplyr::count(df)
  df <- dplyr::ungroup(df)
  # Keep only the maximum per cluster
  df <- dplyr::group_by(df, .data$cluster)
  df <- dplyr::mutate(df, n_max = max(.data$n))
  df <- df[df$n == df$n_max, , drop = FALSE]
  df <- dplyr::ungroup(df)

  # Compute the probability
  df$probability <- df$n / n_samples
  df$n_max <- df$n <- NULL


  # Now collect one set of ranks per cluster
  df <- stats::reshape(
    as.data.frame(df),
    direction = "long",
    varying = setdiff(names(df), c("cluster", "probability")),
    new.row.names = seq_len(prod(dim(df))),
    v.names = "map_ranking",
    timevar = "item",
    times = setdiff(names(df), c("cluster", "probability"))
  )
  df$id <- NULL

  attr(x = df, "reshapeLong") <- NULL # preserves identity to gather() output

  # Sort according to cluster and ranking
  df <- dplyr::arrange(df, .data$cluster, .data$map_ranking)

  if (model_fit$n_clusters[1] == 1) {
    df$cluster <- NULL
  }

  df <- dplyr::as_tibble(df)  # added to solve issue #163. Remove for # 162.

  return(df)

}
