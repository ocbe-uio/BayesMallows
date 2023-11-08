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
#' @param new_rankings Matrix containing the new set of observed rankings of size
#'   n_assessors by n_items.
#' @param n_particles Integer specifying the number of particles.
#' @param augmentation One of "pseudo" and "uniform".
#' @param mcmc_steps Number of Metropolis-Hastings steps to apply in sequential
#'   Monte Carlo.
#' @param ... Optional arguments. Currently not used.
#'
#' @return An updated model, of class "SMCMallows".
#' @export
#'
#'
#' @family modeling
#'
update_mallows <- function(model, new_rankings, ...) {
  UseMethod("update_mallows")
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallows <- function(
    model, new_rankings, n_particles,
    augmentation = "pseudo",
    mcmc_steps = 5, ...) {
  augmentation <- match.arg(augmentation, c("pseudo", "uniform"))
  stopifnot(is.matrix(new_rankings))

  if (augmentation == "pseudo" && !model$metric %in% c("footrule", "spearman")) {
    stop(
      "pseudolikelihood augmentation only possible for metrics 'footrule' ",
      "and 'spearman'."
    )
  }

  alpha_init <- extract_alpha_init(model, n_particles)
  rho_init <- extract_rho_init(model, n_particles)

  ret <- smc_mallows_new_users(
    rankings = t(new_rankings),
    new_rankings = t(new_rankings),
    n_particles = n_particles,
    mcmc_steps = mcmc_steps,
    aug_method = augmentation,
    logz_estimate = model$logz_list$logz_estimate,
    cardinalities = model$logz_list$cardinalities,
    metric = model$metric,
    rho_init = rho_init,
    alpha_init = alpha_init,
    leap_size = floor(model$n_items / 4)
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  ret$alpha <- tidy_parameters$alpha
  ret$rho <- tidy_parameters$rho

  ret$n_items <- model$n_items
  ret$burnin <- 0
  ret$n_clusters <- 1
  ret$rankings <- new_rankings
  ret$n_particles <- n_particles
  ret$mcmc_steps <- mcmc_steps
  ret$logz_list <- model$logz_list
  ret$metric <- model$metric
  ret$items <- model$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}


#' @export
#' @rdname update_mallows
update_mallows.SMCMallows <- function(model, new_rankings, ...) {
  rankings <- rbind(model$rankings, new_rankings)
  alpha_init <- model$alpha_samples
  rho_init <- model$rho_samples
  aug_init <- model$augmented_rankings

  ret <- smc_mallows_new_users(
    rankings = t(rankings),
    new_rankings = t(new_rankings),
    rho_init = rho_init,
    alpha_init = alpha_init,
    n_particles = model$n_particles,
    mcmc_steps = model$mcmc_steps,
    logz_estimate = model$logz_list$logz_estimate,
    cardinalities = model$logz_list$cardinalities,
    metric = model$metric,
    aug_init = aug_init
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  model$alpha <- tidy_parameters$alpha
  model$rho <- tidy_parameters$rho
  model$augmented_rankings <- ret$augmented_rankings
  model$ESS <- ret$ESS
  model$rankings <- rankings

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

