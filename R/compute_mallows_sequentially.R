#' @title Estimate the Bayesian Mallows Model Sequentially
#'
#' @description Compute the posterior distributions of the parameters of the
#' Bayesian Mallows model using sequential Monte Carlo. This is based on the
#' algorithms developed in
#' \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}.
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
#' @return An object of class SMCMallows.
#'
#' @references \insertAllCited{}
#' @export
#'
#' @family modeling
#'
#' @example /inst/examples/compute_mallows_sequentially_examples.R
#'
compute_mallows_sequentially <- function(
    data,
    initial_values,
    model_options = set_model_options(),
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = set_priors(),
    pfun_estimate = NULL) {
  pfun_values <- extract_pfun_values(model_options$metric, data[[1]]$n_items, pfun_estimate)
  validate_class(initial_values, "BayesMallowsPriorSamples")
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
}
