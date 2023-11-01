#' Update a Bayesian Mallows model with new users
#'
#' @param model A model object.
#' @param new_rankings Matrix containing the new set of observed rankings of size
#'   n_assessors by n_items.
#' @param n_particles Integer specifying the number of particles.
#' @param augmentation One of "pseudo" and "uniform".
#' @param mcmc_steps Number of Metropolis-Hastings steps to apply in sequential
#'   Monte Carlo.
#'
#' @return An updated model, of class "SMCMallows".
#' @export
#'
#' @example inst/examples/smc_mallows_new_users_complete_example.R
#'
#' @family modeling
#'
update_mallows <- function(model, new_rankings, n_particles,
                           augmentation = "pseudo",
                           mcmc_steps = 5) {
  UseMethod("update_mallows")
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallows <- function(
    model, new_rankings, n_particles,
    augmentation = "pseudo",
    mcmc_steps = 5) {

  augmentation <- match.arg(augmentation, c("pseudo", "uniform"))
  stopifnot(is.matrix(new_rankings))

  if(any(is.na(new_rankings))) {
    type <- "partial"
  } else {
    type <- "complete"
  }

  alpha_init <- extract_alpha_init(model, n_particles, model$burnin)
  rho_init <- extract_rho_init(model, n_particles, model$burnin)

  ret <- smc_mallows_new_users(
    rankings = t(new_rankings),
    type = type,
    n_particles = n_particles,
    mcmc_steps = mcmc_steps,
    num_new_obs = nrow(new_rankings),
    aug_method = augmentation,
    logz_estimate = model$logz_list$logz_estimate,
    cardinalities = model$logz_list$cardinalities,
    metric = model$metric,
    rho_init = rho_init,
    alpha_init = alpha_init
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  ret$alpha <- tidy_parameters$alpha
  ret$rho <- tidy_parameters$rho
  ret$alpha_samples <- NULL
  ret$rho_samples <- NULL

  ret$n_items <- model$n_items
  ret$burnin <- 0
  ret$n_clusters <- 1
  ret$rankings <- new_rankings
  ret$n_particles <- n_particles
  ret$type <- type
  ret$mcmc_steps <- mcmc_steps
  ret$logz_list <- model$logz_list
  ret$metric <- model$metric
  ret$items <- model$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}


#' @export
#' @rdname update_mallows
update_mallows.SMCMallows <- function(model, new_rankings) {

  rankings <- rbind(model$rankings, new_rankings)
  alpha_init <- extract_alpha_init(model, model$n_particles)
  rho_init <- extract_rho_init(model, model$n_particles)


  ret <- smc_mallows_new_users(
    rankings = t(rankings),
    type = model$type,
    n_particles = model$n_particles,
    mcmc_steps = model$mcmc_steps,
    num_new_obs = nrow(new_rankings),
    logz_estimate = model$logz_list$logz_estimate,
    cardinalities = model$logz_list$cardinalities,
    metric = model$metric,
    rho_init = rho_init,
    alpha_init = alpha_init,
    num_obs = nrow(model$rankings)
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
  result$alpha <- tidy_alpha(
    t(ret$alpha_samples[, ncol(ret$alpha_samples)]), 1, 1)

  rho_mat <- array(dim = c(dim(ret$rho_samples)[[2]], 1, dim(ret$rho_samples)[[1]]))
  rho_mat[, 1, ] <- aperm(ret$rho_samples, c(3, 2, 1))[dim(ret$rho_samples)[[3]], , , drop = TRUE]
  result$rho <- tidy_rho(rho_mat, 1, 1, items)

  result
}

extract_alpha_init <- function(model, n_particles, burnin = 0) {
  sample(
    model$alpha$value[model$alpha$iteration > burnin],
    size = n_particles, replace = TRUE)
}

extract_rho_init <- function(model, n_particles, burnin = 0) {
  rho_init <- reshape(
    data = model$rho[model$rho$iteration > burnin, ],
    v.names = "value",
    idvar = c("chain", "iteration"),
    timevar = "item",
    direction = "wide")
  rho_init <- rho_init[sample(nrow(rho_init), n_particles, replace = TRUE), ]
  rho_init <- rho_init[, grep("^value", colnames(rho_init))]
  t(as.matrix(rho_init))
}
