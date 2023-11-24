#' Update a Bayesian Mallows model with new users
#'
#' Update a Bayesian Mallows model estimated using the Metropolis-Hastings
#' algorithm in [compute_mallows()] using the sequential Monte Carlo
#' algorithm described in
#' \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}. This
#' is useful when new data arrives, and is typically more computationally
#' efficient than running Metropolis-Hastings from scratch.
#'
#' @param model A model object.
#' @param new_data Object returned from [setup_rank_data()] containing new data.
#' @param smc_options SMC specific options returned from [set_smc_options()].
#' @param compute_options An object of class `"BayesMallowsComputeOptions"`
#'   returned from [set_compute_options()].
#' @param priors An object of class `"BayesMallowsPriors"` returned from
#'   [set_priors()]. Defaults to the priors used in `model`.
#' @param ... Optional arguments. Currently not used.
#'
#' @return An updated model, of class "SMCMallows".
#' @export
#'
#'
#' @family modeling
#'
update_mallows <- function(model, new_data, ...) {
  UseMethod("update_mallows")
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallows <- function(
    model, new_data,
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = model$priors,
    ...) {


  if (smc_options$aug_method == "pseudo" && !model$metric %in% c("footrule", "spearman")) {
    stop(
      "pseudolikelihood aug_method only possible for metrics 'footrule' ",
      "and 'spearman'."
    )
  }

  alpha_init <- extract_alpha_init(model, smc_options$n_particles)
  rho_init <- extract_rho_init(model, smc_options$n_particles)

  ret <- smc_mallows_new_users(
    data = new_data,
    new_data = new_data,
    smc_options = smc_options,
    compute_options = compute_options,
    initial_values = list(alpha_init = alpha_init, rho_init = rho_init,
                          aug_init = NULL),
    logz_list = model$logz_list,
    metric = model$metric
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  ret$alpha <- tidy_parameters$alpha
  ret$rho <- tidy_parameters$rho

  ret$smc_options <- smc_options
  ret$compute_options <- compute_options
  ret$n_items <- model$n_items
  ret$burnin <- 0
  ret$n_clusters <- 1
  ret$data <- new_data
  ret$logz_list <- model$logz_list
  ret$metric <- model$metric
  ret$items <- model$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}


#' @export
#' @rdname update_mallows
update_mallows.SMCMallows <- function(model, new_data, ...) {

  data <- setup_rank_data(rankings = rbind(model$data$rankings, new_data$rankings))
  alpha_init <- model$alpha_samples
  rho_init <- model$rho_samples
  aug_init <- model$augmented_rankings

  ret <- smc_mallows_new_users(
    data = data,
    new_data = new_data,
    smc_options = model$smc_options,
    compute_options = model$compute_options,
    initial_values = list(alpha_init = alpha_init, rho_init = rho_init,
                          aug_init = aug_init),
    logz_list = model$logz_list,
    metric = model$metric
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  model$alpha <- tidy_parameters$alpha
  model$rho <- tidy_parameters$rho
  model$augmented_rankings <- ret$augmented_rankings
  model$ESS <- ret$ESS
  model$data <- data

  class(model) <- c("SMCMallows", "BayesMallows")
  model
}


tidy_smc <- function(ret, items) {
  result <- list()
  result$alpha <- tidy_alpha(matrix(ret$alpha_samples, nrow = 1), 1, 1)

  rho_mat <- array(dim = c(dim(ret$rho_samples)[[1]], 1, dim(ret$rho_samples)[[2]]))
  rho_mat[, 1, ] <- ret$rho_samples
  result$rho <- tidy_rho(rho_mat, 1, 1, items)

  result
}

extract_alpha_init <- function(model, n_particles) {
  thinned_inds <- floor(
    seq(from = model$burnin + 1, to = ncol(model$alpha_samples),
        length.out = n_particles))
  model$alpha_samples[1,  thinned_inds, drop = TRUE]
}

extract_rho_init <- function(model, n_particles) {
  thinned_inds <- floor(
    seq(from = model$burnin + 1, to = dim(model$rho_samples)[[3]],
        length.out = n_particles))
  model$rho_samples[, 1, thinned_inds, drop = TRUE]
}

