#' @title Estimate the Bayesian Mallows Model Sequentially
#'
#' @description Compute the posterior distributions of the parameters of the
#'   Bayesian Mallows model using sequential Monte Carlo. This is based on the
#'   algorithms developed in
#'   \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}.
#'   This function differs from [update_mallows()] in that it takes all the data
#'   at once, and uses SMC to fit the model step-by-step. Used in this way, SMC
#'   is an alternative to Metropolis-Hastings, which may work better in some
#'   settings. In addition, it allows visualization of the learning process.
#'
#' @param data A list of objects of class "BayesMallowsData" returned from
#'   [setup_rank_data()]. Each list element is interpreted as the data belonging
#'   to a given timepoint.
#' @param initial_values An object of class "BayesMallowsPriorSamples" returned
#'   from [sample_prior()].
#' @param model_options An object of class "BayesMallowsModelOptions" returned
#'   from [set_model_options()].
#' @param smc_options An object of class "SMCOptions" returned from
#'   [set_smc_options()].
#' @param compute_options An object of class "BayesMallowsComputeOptions"
#'   returned from [set_compute_options()].
#' @param priors An object of class "BayesMallowsPriors" returned from
#'   [set_priors()].
#'
#' @param pfun_estimate Object returned from [estimate_partition_function()].
#'   Defaults to \code{NULL}, and will only be used for footrule, Spearman, or
#'   Ulam distances when the cardinalities are not available, cf.
#'   [get_cardinalities()].
#'
#' @return An object of class BayesMallowsSequential.
#'
#' @details This function is very new, and plotting functions and other tools
#' for visualizing the posterior distribution do not yet work. See the examples
#' for some workarounds.
#'
#'
#' @references \insertAllCited{}
#' @export
#'
#' @family modeling
#'
#' @example /inst/examples/compute_mallows_sequentially_example.R
#'
compute_mallows_sequentially <- function(
    data,
    initial_values,
    model_options = set_model_options(),
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = set_priors(),
    pfun_estimate = NULL) {
  validate_class(initial_values, "BayesMallowsPriorSamples")
  if (!is.list(data) | !all(vapply(data, inherits, logical(1), "BayesMallowsData"))) {
    stop("data must be a list of BayesMallowsData objects.")
  }
  if (any(
    vapply(data, function(x) {
      is.null(x$user_ids) || length(x$user_ids) == 0
    }, logical(1))
  )) {
    stop("User IDs must be set.")
  }
  data <- lapply(data, function(x) {
    if(!is.null(x$preferences)) {
      x$preferences <- as.matrix(x$preferences)
    } else {
      x$preferences <- matrix(0, 0, 0)
    }
    x
  })
  pfun_values <- extract_pfun_values(model_options$metric, data[[1]]$n_items, pfun_estimate)
  alpha_init <- sample(initial_values$alpha, smc_options$n_particles, replace = TRUE)
  rho_init <- initial_values$rho[, sample(ncol(initial_values$rho), smc_options$n_particles, replace = TRUE)]

  ret <- run_smc(
    data = flush(data[[1]]),
    new_data = data,
    model_options = model_options,
    smc_options = smc_options,
    compute_options = compute_options,
    priors = priors,
    initial_values = list(alpha_init = alpha_init, rho_init = rho_init, aug_init = NULL),
    pfun_values = pfun_values,
    pfun_estimate = pfun_estimate
  )
  class(ret) <- "SMCMallows"
  ret
}
