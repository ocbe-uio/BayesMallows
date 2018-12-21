#' Estimate Partition Function
#'
#' Estimate the logarithm of the partition function of the Mallows rank model. Choose
#' between the importance sampling algorithm described in
#' \insertCite{vitelli2018}{BayesMallows} and the IPFP algorithm for computing an
#' asymptotic approximation described in \insertCite{mukherjee2016}{BayesMallows}.
#'
#' @param method Character string specifying the method to use in order to
#'   estimate the logarithm of the partition function. Available options are
#'   \code{"importance_sampling"} and \code{"asymptotic"}.
#'
#' @param alpha_vector Numeric vector of \eqn{\alpha} values over which to
#'   compute the importance sampling estimate.
#'
#' @param n_items Integer specifying the number of items.
#'
#' @param metric Character string specifying the distance measure to use.
#'   Available options are \code{"footrule"} and \code{"spearman"}
#'
#' @param nmc Integer specifying the number of Monte Carlo samples to use in the
#'   importance sampling. Only used when \code{method = "importance_sampling"}.
#'
#' @param degree Integer specifying the degree of the polynomial used to
#'   estimate \eqn{\log(\alpha)} from the grid of values provided by the
#'   importance sampling estimate.
#'
#' @param n_iterations Integer specifying the number of iterations to use in the
#'   asymptotic approximation of the partition function. Only used when
#'   \code{method = "asymptotic"}.
#'
#' @param K Integer specifying the parameter \eqn{K} in the
#' asymptotic approximation of the partition function. Only used when
#' \code{method = "asymptotic"}.
#'
#' @return A vector of length \code{degree} which can be supplied to the
#'   \code{logz_estimate} argument of \code{\link{compute_mallows}}.
#'
#' @param cl Optional computing cluster used for parallelization, returned
#' from \code{parallel::makeCluster}. Defaults to \code{NULL}. Only used when
#' \code{method = "importance_sampling"}.
#'
#' @export
#'
#' @references \insertAllCited{}
#'
#' @example /inst/examples/estimate_partition_function_example.R
#'
estimate_partition_function <- function(method = "importance_sampling",
                                        alpha_vector, n_items, metric,
                                        nmc, degree, n_iterations, K, cl = NULL){

  stopifnot(degree < length(alpha_vector))

  if(method == "importance_sampling"){
    if(!is.null(cl)){
      # Split nmc into each cluster
      nmc_vec <- rep(floor(nmc / length(cl)), length(cl))
      i <- 1
      while(sum(nmc_vec) != nmc){
        nmc_vec[i] <- nmc_vec[i] + 1
        if(i > length(cl)) break()
      }
      parallel::clusterExport(cl, c("alpha_vector", "n_items", "metric"),
                              envir = environment())
      estimates <- parallel::parLapply(cl, nmc_vec, function(x){
        compute_importance_sampling_estimate(alpha_vector = alpha_vector, n_items = n_items,
                                             metric = metric, nmc = x)
      })
      log_z <- purrr::map2_dbl(estimates[[1]], estimates[[2]], mean)
    } else {
      log_z <- as.numeric(
        compute_importance_sampling_estimate(
          alpha_vector = alpha_vector, n_items = n_items,
          metric = metric, nmc = nmc))
    }

    # Compute the estimate at each discrete alpha value
    estimate <- dplyr::tibble(
      alpha = alpha_vector,
      log_z = log_z
    )
  } else if(method == "asymptotic"){
    estimate <- dplyr::tibble(
      alpha = alpha_vector,
      log_z = as.numeric(
        asymptotic_partition_function(
          alpha_vector = alpha_vector, n_items = n_items,
          metric = metric, K = K, n_iterations = n_iterations))
    )
  }

  # Fit a regression model
  form <- stats::as.formula(paste("log_z ~ ",
                           paste("I( alpha^", seq(from = 1, to = degree, by = 1), ")",
                                 collapse = "+")))
  stats::lm(form, data = estimate)$coefficients

}
