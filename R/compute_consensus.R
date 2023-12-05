#' @title Compute Consensus Ranking
#' @description Compute the consensus ranking using either cumulative
#'   probability (CP) or maximum a posteriori (MAP) consensus
#'   \insertCite{vitelli2018}{BayesMallows}. For mixture models, the consensus
#'   is given for each mixture. Consensus of augmented ranks can also be
#'   computed for each assessor, by setting `parameter = "Rtilde"`.
#'
#' @param model_fit A model fit.
#' @param type Character string specifying which consensus to compute. Either
#'   `"CP"` or `"MAP"`. Defaults to `"CP"`.
#' @param burnin A numeric value specifying the number of iterations to discard
#'   as burn-in. Defaults to `model_fit$burnin`, and must be provided if
#'   `model_fit$burnin` does not exist. See [assess_convergence()].
#' @param parameter Character string defining the parameter for which to compute
#'   the consensus. Defaults to `"rho"`. Available options are `"rho"` and
#'   `"Rtilde"`, with the latter giving consensus rankings for augmented ranks.
#' @param assessors When `parameter = "rho"`, this integer vector is used to
#'   define the assessors for which to compute the augmented ranking. Defaults
#'   to `1L`, which yields augmented rankings for assessor 1.
#' @param ... Other arguments passed on to other methods. Currently not used.
#'
#' @references \insertAllCited{}
#' @export
#' @example /inst/examples/compute_consensus_example.R
#'
#' @family posterior quantities
#'
compute_consensus <- function(model_fit, ...) {
  UseMethod("compute_consensus")
}

#' @export
#' @rdname compute_consensus
compute_consensus.BayesMallows <- function(
    model_fit, type = c("CP", "MAP"), burnin = model_fit$burnin,
    parameter = c("rho", "Rtilde"), assessors = 1L, ...) {
  if (is.null(burnin)) stop("Please specify the burnin.")
  stopifnot(burnin < model_fit$nmc)
  type <- match.arg(type, c("CP", "MAP"))
  parameter <- match.arg(parameter, c("rho", "Rtilde"))

  if (parameter == "Rtilde" &&
    !inherits(model_fit$augmented_data, "data.frame")) {
    stop("For augmented ranks, please refit model with option 'save_aug = TRUE'.")
  }

  if (parameter == "rho") {
    df <- model_fit$rho[model_fit$rho$iteration > burnin, , drop = FALSE]
    if (type == "CP") {
      df <- cpc_bm(df)
    } else if (type == "MAP") {
      df <- cpm_bm(df)
    }
  } else if (parameter == "Rtilde") {
    df <- model_fit$augmented_data[
      model_fit$augmented_data$iteration > burnin &
      model_fit$augmented_data$assessor %in% assessors, , drop = FALSE]

    names(df)[names(df) == "assessor"] <- "cluster"
    class(df) <- c("consensus_BayesMallows", "tbl_df", "tbl", "data.frame")

    df <- if (type == "CP") {
      df <- cpc_bm(df)
    } else if (type == "MAP") {
      df <- cpm_bm(df)
    }

    if ("cluster" %in% names(df)) {
      df$assessor <- as.numeric(df$cluster)
      df$cluster <- NULL
    }
  }

  row.names(df) <- NULL
  as.data.frame(df)
}


#' @export
#' @rdname compute_consensus
compute_consensus.SMCMallows <- function(
    model_fit, type = c("CP", "MAP"), parameter = "rho", ...) {
  parameter <- match.arg(parameter, "rho")
  model_fit$burnin <- 0
  model_fit$nmc <- model_fit$n_particles
  NextMethod("compute_consensus")
}


# Internal function for finding CP consensus.
find_cpc <- function(group_df, group_var = "cluster") {
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

    if (nrow(tmp_df) >= 1) {
      # Keep the max only. This filtering must be done after the first filter,
      # since we take the maximum among the filtered values
      tmp_df <- do.call(
        rbind,
        lapply(split(tmp_df, f = tmp_df[group_var]), function(x) {
          x[x$cumprob == max(x$cumprob), ]
        })
      )
      # Add the ranking
      tmp_df$ranking <- i

      # Select the columns we want to keep, and put them in result
      result <- rbind(
        result,
        tmp_df[, c("cluster", "ranking", "item", "cumprob"), drop = FALSE]
      )
    }
  }
  return(result)
}

aggregate_cp_consensus <- function(df) {
  # Convert items and cluster to character, since factor levels are not needed in this case
  df$item <- as.character(df$item)
  df$cluster <- as.character(df$cluster)

  df <- aggregate(
    list(n = df$iteration),
    by = list(
      item = as.character(df$item),
      cluster = as.character(df$cluster), value = df$value
    ),
    FUN = length
  )

  # Arrange according to value, per item and cluster
  do.call(rbind, lapply(split(df, f = ~ item + cluster), function(x) {
    x <- x[order(x$value), ]
    x$cumprob <- cumsum(x$n) / sum(x$n)
    x
  }))
}

aggregate_map_consensus <- function(df, n_samples) {
  # Group by everything except iteration, and count the unique combinations
  df <- aggregate(list(n = df$iteration), df[, setdiff(names(df), "iteration")],
    FUN = length
  )
  # Keep only the maximum per cluster
  df <- do.call(rbind, lapply(split(df, f = df$cluster), function(x) {
    x$n_max <- max(x$n)
    x[x$n == x$n_max, , drop = FALSE]
  }))

  # Compute the probability
  df$probability <- df$n / n_samples
  df$n_max <- df$n <- NULL
  df
}

cpc_bm <- function(df) {
  df <- aggregate_cp_consensus(df)
  df <- find_cpc(df)
  df[order(df$cluster, df$ranking), ]
}

cpm_bm <- function(df) {
  n_samples <- length(unique(df$iteration))

  # Reshape to get items along columns
  df <- stats::reshape(as.data.frame(df),
    direction = "wide",
    idvar = c("chain", "cluster", "iteration"),
    timevar = "item"
  )
  df$chain <- NULL
  names(df) <- gsub("^value\\.", "", names(df))

  df <- aggregate_map_consensus(df, n_samples)

  # Now collect one set of ranks per cluster
  df$id <- seq_len(nrow(df))
  df <- stats::reshape(as.data.frame(df),
    direction = "long",
    varying = setdiff(names(df), c("cluster", "probability", "id")),
    v.names = "map_ranking",
    timevar = "item",
    idvar = c("cluster", "probability", "id"),
    times = setdiff(names(df), c("cluster", "probability", "id"))
  )
  rownames(df) <- NULL
  df$id <- NULL

  # Sort according to cluster and ranking
  df[order(df$cluster, df$map_ranking),
    c("cluster", "map_ranking", "item", "probability"),
    drop = FALSE
  ]
}
