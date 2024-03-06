#' @title Setup rank data
#'
#' @description Prepare rank or preference data for further analyses.
#'
#' @param rankings A matrix of ranked items, of size `n_assessors x n_items`.
#'   See [create_ranking()] if you have an ordered set of items that need to be
#'   converted to rankings. If `preferences` is provided, `rankings` is an
#'   optional initial value of the rankings. If `rankings` has column names,
#'   these are assumed to be the names of the items. `NA` values in rankings are
#'   treated as missing data and automatically augmented; to change this
#'   behavior, see the `na_action` argument to [set_model_options()]. A vector
#'   length `n_items` is silently converted to a matrix of length `1 x n_items`,
#'   and names (if any), are used as column names.
#'
#' @param preferences A data frame with one row per pairwise comparison, and
#'   columns `assessor`, `top_item`, and `bottom_item`. Each column contains the
#'   following:
#' \itemize{
#' \item `assessor` is a numeric vector containing the assessor index.
#'
#' \item `bottom_item` is a numeric vector containing the index of the item that
#'       was disfavored in each pairwise comparison.
#'
#' \item `top_item` is a numeric vector containing the index of the item that was
#'       preferred in each pairwise comparison.
#' }
#'   So if we have two assessors and five items, and assessor 1 prefers item 1
#'   to item 2 and item 1 to item 5, while assessor 2 prefers item 3 to item 5,
#'   we have the following `df`:
#' \tabular{rrr}{
#' **assessor** \tab **bottom_item** \tab **top_item**\cr
#' 1 \tab 2 \tab 1\cr
#' 1 \tab 5 \tab 1\cr
#' 2 \tab 5 \tab 3\cr
#' }
#'
#' @param user_ids Optional `numeric` vector of user IDs. Only only used by
#'   [update_mallows()]. If provided, new data can consist of updated partial
#'   rankings from users already in the dataset, as described in Section 6 of
#'   \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}.
#'
#' @param observation_frequency A vector of observation frequencies (weights) to
#'   apply do each row in `rankings`. This can speed up computation if a large
#'   number of assessors share the same rank pattern. Defaults to `NULL`, which
#'   means that each row of `rankings` is multiplied by 1. If provided,
#'   `observation_frequency` must have the same number of elements as there are
#'   rows in `rankings`, and `rankings` cannot be `NULL`. See
#'   [compute_observation_frequency()] for a convenience function for computing
#'   it.
#'
#' @param validate_rankings Logical specifying whether the rankings provided (or
#'   generated from `preferences`) should be validated. Defaults to `TRUE`.
#'   Turning off this check will reduce computing time with a large number of
#'   items or assessors.
#'
#' @param na_action Character specifying how to deal with `NA` values in the
#'   `rankings` matrix, if provided. Defaults to `"augment"`, which means that
#'   missing values are automatically filled in using the Bayesian data
#'   augmentation scheme described in
#'   \insertCite{vitelli2018;textual}{BayesMallows}. The other options for this
#'   argument are `"fail"`, which means that an error message is printed and the
#'   algorithm stops if there are `NA`s in `rankings`, and `"omit"` which simply
#'   deletes rows with `NA`s in them.
#'
#' @param cl Optional computing cluster used for parallelization when generating
#'   transitive closure based on preferences, returned from
#'   [parallel::makeCluster()]. Defaults to `NULL`.
#'
#' @param shuffle_unranked Logical specifying whether or not to randomly permute
#'   unranked items in the initial ranking. When `shuffle_unranked=TRUE` and
#'   `random=FALSE`, all unranked items for each assessor are randomly permuted.
#'   Otherwise, the first ordering returned by `igraph::topo_sort()` is
#'   returned.
#'
#' @param random Logical specifying whether or not to use a random initial
#'   ranking. Defaults to `FALSE`. Setting this to `TRUE` means that all
#'   possible orderings consistent with the stated pairwise preferences are
#'   generated for each assessor, and one of them is picked at random.
#'
#' @param random_limit Integer specifying the maximum number of items allowed
#'   when all possible orderings are computed, i.e., when `random=TRUE`.
#'   Defaults to `8L`.
#'
#' @param timepoint Integer vector specifying the timepoint. Defaults to `NULL`,
#'   which means that a vector of ones, one for each observation, is generated.
#'   Used by [update_mallows()] to identify data with a given iteration of the
#'   sequential Monte Carlo algorithm. If not `NULL`, must contain one integer
#'   for each row in `rankings`.
#'
#' @param n_items Integer specifying the number of items. Defaults to `NULL`,
#'   which means that the number of items is inferred from `rankings` or from
#'   `preferences`. Setting `n_items` manually can be useful with pairwise
#'   preference data in the SMC algorithm, i.e., when `rankings` is `NULL` and
#'   `preferences` is non-`NULL`, and contains a small number of pairwise
#'   preferences for a subset of users and items.
#'
#' @note Setting `random=TRUE` means that all possible orderings of each
#'   assessor's preferences are generated, and one of them is picked at random.
#'   This can be useful when experiencing convergence issues, e.g., if the MCMC
#'   algorithm does not mix properly. However, finding all possible orderings is
#'   a combinatorial problem, which may be computationally very hard. The result
#'   may not even be possible to fit in memory, which may cause the R session to
#'   crash. When using this option, please try to increase the size of the
#'   problem incrementally, by starting with smaller subsets of the complete
#'   data. An example is given below.
#'
#'   It is assumed that the items are labeled starting from 1. For example, if a
#'   single comparison of the following form is provided, it is assumed that
#'   there is a total of 30 items (`n_items=30`), and the initial ranking is a
#'   permutation of these 30 items consistent with the preference 29<30.
#'
#' \tabular{rrr}{
#' **assessor** \tab **bottom_item** \tab **top_item**\cr
#' 1 \tab 29 \tab 30\cr
#' }
#'
#'   If in reality there are only two items, they should be relabeled to 1 and
#'   2, as follows:
#'
#' \tabular{rrr}{
#' **assessor** \tab **bottom_item** \tab **top_item**\cr
#' 1 \tab 1 \tab 2\cr
#' }
#'
#'
#'
#' @return An object of class `"BayesMallowsData"`, to be provided in the `data`
#'   argument to [compute_mallows()].
#' @export
#'
#' @family preprocessing
#'
#' @references \insertAllCited{}
#'
setup_rank_data <- function(
    rankings = NULL,
    preferences = NULL,
    user_ids = numeric(),
    observation_frequency = NULL,
    validate_rankings = TRUE,
    na_action = c("augment", "fail", "omit"),
    cl = NULL,
    shuffle_unranked = FALSE,
    random = FALSE,
    random_limit = 8L,
    timepoint = NULL,
    n_items = NULL) {
  na_action <- match.arg(na_action, c("augment", "fail", "omit"))
  if (!is.null(rankings) && !is.null(n_items)) {
    stop("n_items can only be set when rankings=NULL")
  }

  if (is.null(rankings) && is.null(preferences)) {
    stop("Either rankings or preferences (or both) must be provided.")
  }

  preferences <- generate_transitive_closure(preferences, cl)
  if (!is.null(rankings)) {
    if (na_action == "fail" && any(is.na(rankings))) {
      stop("rankings matrix contains NA values")
    }
    if (!is.matrix(rankings)) {
      rankings <- matrix(rankings,
        nrow = 1,
        dimnames = list(NULL, names(rankings))
      )
    }

    if (na_action == "omit" && any(is.na(rankings))) {
      keeps <- apply(rankings, 1, function(x) !any(is.na(x)))
      print(paste(
        "Omitting", sum(!keeps),
        "row(s) from rankings due to NA values"
      ))
      rankings <- rankings[keeps, , drop = FALSE]
    }
  } else {
    if (is.null(n_items)) n_items <- max(preferences[, c("bottom_item", "top_item")])
    rankings <- generate_initial_ranking(
      preferences, n_items, cl, shuffle_unranked, random, random_limit
    )
  }

  if (!is.null(observation_frequency)) {
    validate_positive_vector(observation_frequency)
    if (nrow(rankings) != length(observation_frequency)) {
      stop(
        "observation_frequency must be of same ",
        "length as the number of rows in rankings"
      )
    }
  } else {
    observation_frequency <- rep(1, nrow(rankings))
  }

  if (validate_rankings &&
    !all(apply(rankings, 1, validate_permutation))) {
    stop("invalid permutations provided in rankings matrix")
  }
  n_items <- ncol(rankings)

  if (!is.null(colnames(rankings))) {
    items <- colnames(rankings)
  } else {
    items <- paste("Item", seq(from = 1, to = n_items, by = 1))
  }

  if (is.null(timepoint)) timepoint <- rep(1, nrow(rankings))
  if (length(timepoint) != nrow(rankings)) {
    stop("must have one timepoint per row")
  }

  constraints <- generate_constraints(preferences, n_items, cl)
  consistent <- matrix(integer(0))
  n_assessors <- nrow(rankings)
  any_missing <- any(is.na(rankings))
  augpair <- !is.null(preferences)
  stopifnot(is.numeric(user_ids))

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsData"
  ret
}
