#' Estimate Partition Function
#'
#' Estimate the logarithm of the partition function of the Mallows rank model.
#' Choose between the importance sampling algorithm described in
#' \insertCite{vitelli2018}{BayesMallows} and the IPFP algorithm for computing
#' an asymptotic approximation described in
#' \insertCite{mukherjee2016}{BayesMallows}.
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
#' @return A vector of length `degree` which can be supplied to the
#'   `logz_estimate` argument of [compute_mallows()].
#'
#' @param cl Optional computing cluster used for parallelization, returned from
#'   [parallel::makeCluster()]. Defaults to `NULL`. Only used when `method =
#'   "importance_sampling"`.
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
      # Split n_iterations into each cluster
      n_iterations_vec <- rep(floor(n_iterations / length(cl)), length(cl))
      i <- 1
      while (sum(n_iterations_vec) != n_iterations) {
        n_iterations_vec[i] <- n_iterations_vec[i] + 1
        if (i > length(cl)) break
      }
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
    estimate <- data.frame(
      alpha = alpha_vector,
      log_z = log_z
    )
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

  form <- stats::as.formula(paste(
    "log_z ~ ",
    paste("I( alpha^", seq(from = 1, to = degree, by = 1), ")",
      collapse = "+"
    )
  ))
  stats::lm(form, data = estimate)$coefficients
}

prepare_partition_function <- function(logz_estimate = NULL, metric, n_items) {
  ret <- list(cardinalities = NULL, logz_estimate = NULL)
  class(ret) <- "BayesMallowsPartitionFunction"

  if (!is.null(logz_estimate)) {
    ret$logz_estimate <- logz_estimate
    return(ret)
  }

  if (metric %in% c("cayley", "hamming", "kendall")) {
    return(ret)
  }

  tryCatch(
    {
      relevant_params <- get_cardinalities(n_items, metric)
      ret$cardinalities <- relevant_params$value
      return(ret)
    },
    error = function(e) {
      stop("Partition function not available. Please compute an estimate using estimate_partition_function().")
    }
  )
}
