#' @title Estimate Partition Function
#'
#' @description
#' Estimate the logarithm of the partition function of the Mallows rank model.
#' Choose between the importance sampling algorithm described in
#' \insertCite{vitelli2018}{BayesMallows} and the IPFP algorithm for computing
#' an asymptotic approximation described in
#' \insertCite{mukherjee2016}{BayesMallows}. Note that exact partition functions
#' can be computed efficiently for Cayley, Hamming and Kendall distances with
#' any number of items, for footrule distances with up to 50 items, Spearman
#' distance with up to 20 items, and Ulam distance with up to 60 items. This
#' function is thus intended for the complement of these cases. See
#' [get_cardinalities()] for details.
#'
#' @param method Character string specifying the method to use in order to
#'   estimate the logarithm of the partition function. Available options are
#'   `"importance_sampling"` and `"asymptotic"`.
#'
#' @param alpha_vector Numeric vector of \eqn{\alpha} values over which to
#'   compute the importance sampling estimate.
#'
#' @param n_items Integer specifying the number of items.
#'
#' @param metric Character string specifying the distance measure to use.
#'   Available options are `"footrule"` and `"spearman"` when `method =
#'   "asymptotic"` and in addition `"cayley"`, `"hamming"`, `"kendall"`, and
#'   `"ulam"` when `method = "importance_sampling"`.
#'
#' @param n_iterations Integer specifying the number of iterations to use. When
#'   `method = "importance_sampling"`, this is the number of Monte Carlo samples
#'   to generate. When `method = "asymptotic"`, on the other hand, it represents
#'   the number of iterations of the IPFP algorithm.
#'
#' @param K Integer specifying the parameter \eqn{K} in the asymptotic
#'   approximation of the partition function. Only used when `method =
#'   "asymptotic"`. Defaults to 20.
#'
#' @param cl Optional computing cluster used for parallelization, returned from
#'   [parallel::makeCluster()]. Defaults to `NULL`. Only used when `method =
#'   "importance_sampling"`.
#'
#'
#' @return A matrix with two column and number of rows equal the degree of the
#'   fitted polynomial approximating the partition function. The matrix can be
#'   supplied to the `pfun_estimate` argument of [compute_mallows()].
#'
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/estimate_partition_function_example.R
#' @family partition function
#'
estimate_partition_function <- function(
    method = c("importance_sampling", "asymptotic"),
    alpha_vector, n_items, metric,
    n_iterations, K = 20, cl = NULL) {
  degree <- min(10, length(alpha_vector))

  method <- match.arg(method, c("importance_sampling", "asymptotic"))

  if (method == "importance_sampling") {
    metric <- match.arg(metric, c(
      "footrule", "spearman", "cayley", "hamming",
      "kendall", "ulam"
    ))
    if (!is.null(cl)) {
      n_iterations_vec <- count_jobs_per_cluster(n_iterations, length(cl))
      parallel::clusterExport(cl, c("alpha_vector", "n_items", "metric"),
        envir = environment()
      )
      parallel::clusterSetRNGStream(cl)

      estimates <- parallel::parLapply(cl, n_iterations_vec, function(x) {
        compute_importance_sampling_estimate(
          alpha_vector = alpha_vector, n_items = n_items,
          metric = metric, nmc = x
        )
      })

      log_z <- rowMeans(do.call(cbind, estimates))
    } else {
      log_z <- as.numeric(
        compute_importance_sampling_estimate(
          alpha_vector = alpha_vector, n_items = n_items,
          metric = metric, nmc = n_iterations
        )
      )
    }

    # Compute the estimate at each discrete alpha value
    estimate <- data.frame(alpha = alpha_vector, log_z = log_z)
  } else if (method == "asymptotic") {
    metric <- match.arg(metric, c("footrule", "spearman"))

    estimate <- data.frame(
      alpha = alpha_vector,
      log_z = as.numeric(
        asymptotic_partition_function(
          alpha_vector = alpha_vector, n_items = n_items,
          metric = metric, K = K, n_iterations = n_iterations
        )
      )
    )
  }

  power <- seq(from = 0, to = degree, by = 1)
  form <- stats::as.formula(paste(
    "log_z ~ 0 + ", paste("I( alpha^", power, ")", collapse = "+")
  ))
  matrix(c(power, stats::lm(form, data = estimate)$coefficients), ncol = 2)
}

prepare_partition_function <- function(metric, n_items) {
  if (metric %in% c("cayley", "hamming", "kendall")) {
    return(NULL)
  }

  tryCatch(
    return(as.matrix(get_cardinalities(n_items, metric))),
    error = function(e) {
      stop(
        "Partition function not available. ",
        "Please compute an estimate using estimate_partition_function()."
      )
    }
  )
}
