#' Update a Bayesian Mallows model with new users
#'
#' Update a Bayesian Mallows model estimated using the Metropolis-Hastings
#' algorithm in [compute_mallows()] using the sequential Monte Carlo algorithm
#' described in
#' \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}.
#'
#' @param model A model object of class "BayesMallows" returned from
#'   [compute_mallows()], an object of class "SMCMallows" returned from this
#'   function, or an object of class "BayesMallowsPriorSamples" returned from
#'   [sample_prior()].
#' @param new_data An object of class "BayesMallowsData" returned from
#'   [setup_rank_data()]. The object should contain the new data being provided.
#' @param model_options An object of class "BayesMallowsModelOptions" returned
#'   from [set_model_options()].
#' @param smc_options An object of class "SMCOptions" returned from
#'   [set_smc_options()].
#' @param compute_options An object of class "BayesMallowsComputeOptions"
#'   returned from [set_compute_options()].
#' @param priors An object of class "BayesMallowsPriors" returned from
#'   [set_priors()]. Defaults to the priors used in `model`.
#' @param pfun_estimate Object returned from [estimate_partition_function()].
#'   Defaults to \code{NULL}, and will only be used for footrule, Spearman, or
#'   Ulam distances when the cardinalities are not available, cf.
#'   [get_cardinalities()]. Only used by the specialization for objects of type
#'   "BayesMallowsPriorSamples".
#' @param ... Optional arguments. Currently not used.
#'
#' @return An updated model, of class "SMCMallows".
#' @export
#'
#' @family modeling
#'
#' @example /inst/examples/update_mallows_example.R
#'
update_mallows <- function(model, new_data, ...) {
  UseMethod("update_mallows")
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallowsPriorSamples <- function(
    model, new_data, model_options = set_model_options(),
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = model$priors,
    pfun_estimate = NULL,
    ...) {
  alpha_init <- sample(model$alpha, smc_options$n_particles, replace = TRUE)
  rho_init <- model$rho[, sample(ncol(model$rho), smc_options$n_particles, replace = TRUE) ]
  pfun_values <- extract_pfun_values(model_options, new_data, pfun_estimate)

  run_common_part(
    data = new_data, new_data = new_data, model_options = model_options,
    smc_options = smc_options, compute_options = compute_options,
    priors = priors,
    initial_values = list(
      alpha_init = alpha_init, rho_init = rho_init,
      aug_init = NULL
    ),
    pfun_list = list(pfun_values = pfun_values,
                     pfun_estimate = pfun_estimate),
    model = model)
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallows <- function(
    model, new_data,
    model_options = set_model_options(),
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = model$priors,
    ...) {
  if (is.null(model$burnin)) stop("Burnin must be set.")

  alpha_init <- extract_alpha_init(model, smc_options$n_particles)
  rho_init <- extract_rho_init(model, smc_options$n_particles)

  run_common_part(
    data = new_data, new_data = new_data, model_options = model_options,
    smc_options = smc_options, compute_options = compute_options,
    priors = priors,
    initial_values = list(
      alpha_init = alpha_init, rho_init = rho_init,
      aug_init = NULL
    ),
    pfun_list = list(pfun_values = model$pfun_values,
                     pfun_estimate = model$pfun_estimate),
    model = model)
}

#' @export
#' @rdname update_mallows
update_mallows.SMCMallows <- function(model, new_data, ...) {
  datlist <- prepare_new_data(model, new_data)

  ret <- run_smc(
    data = datlist$data,
    new_data = datlist$new_data,
    model_options = model$model_options,
    smc_options = model$smc_options,
    compute_options = model$compute_options,
    priors = model$priors,
    initial_values = list(
      alpha_init = model$alpha_samples,
      rho_init = model$rho_samples,
      aug_init = model$augmented_rankings
    ),
    pfun_values = model$pfun_values,
    pfun_estimate = model$pfun_estimate
  )
  tidy_parameters <- tidy_smc(ret, model$items)
  model$alpha <- tidy_parameters$alpha
  model$rho <- tidy_parameters$rho
  model$augmented_rankings <- ret$augmented_rankings
  model$data <- datlist$data

  class(model) <- c("SMCMallows", "BayesMallows")
  model
}
