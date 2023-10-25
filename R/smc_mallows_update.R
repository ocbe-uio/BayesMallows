#' Generic function for updating SMC Mallows models
#'
#' @param model A model object.
#' @param rankings A matrix with rankings.
#' @param ... Other optional arguments.
#'
#' @return An updated model.
#' @export
#'
smc_mallows_update <- function(model, rankings, ...) {
  UseMethod("smc_mallows_update")
}

#' Update an SMC Mallows model with new users
#'
#' @param model Object of class \code{c("SMCMallowsNewUsers","SMCMallows")}
#'   returned from \code{\link{smc_mallows_new_users}}.
#' @param rankings Matrix containing the new set of observed rankings of size
#'   n_assessors by n_items.
#' @param timesteps Integer specifying the number of timesteps. Defaults to 1.
#' @param n_particles Integer specifying the number of particles. Defaults to
#'   the number of particles used when computing \code{model}.
#' @param num_new_obs Integer specifying the number of new observations.
#'   Defaults to \code{nrow(rankings)}.
#' @param verbose Logical specifying whether to print out the progress of the
#'   SMC-Mallows algorithm. Defaults to \code{FALSE}.
#' @param ... Optional additional arguments. Currently not used.
#'
#' @return An updated model, of class
#'   \code{c("SMCMallowsNewUsers","SMCMallows")}.
#' @export
#'
#' @example inst/examples/smc_mallows_new_users_complete_example.R
#'
#' @family modeling
#'
smc_mallows_update.SMCMallowsNewUsers <- function(
    model, rankings,
    timesteps = 1,
    num_new_obs = nrow(rankings),
    n_particles = model$n_particles,
    verbose = FALSE,
    ...
    ) {

  if(model$type != "complete") {
    stop("smc_mallows_update currently only implemented for complete data")
  }

  stopifnot(is.matrix(rankings))
  num_new_obs <- num_new_obs # fighting lazy evalution
  n_users <- model$num_obs + nrow(rankings)
  rankings <- rbind(rankings, model$rankings)

  ret <- smc_mallows_new_users_cpp(
    rankings = rankings,
    type = model$type,
    n_items = model$n_items,
    n_users = n_users,
    n_particles = model$n_particles,
    timesteps = timesteps,
    mcmc_kernel_app = model$mcmc_kernel_app,
    num_new_obs = num_new_obs,
    alpha_prop_sd = model$alpha_prop_sd,
    lambda = model$lambda,
    alpha_max = model$alpha_max,
    alpha = model$alpha,
    aug_method = model$aug_method,
    logz_estimate = model$logz_list$logz_estimate,
    cardinalities = model$logz_list$cardinalities,
    verbose = verbose,
    metric = model$metric,
    leap_size = model$leap_size,
    rho_init = model$rho_samples[, , dim(model$rho_samples)[[3]]],
    alpha_init = model$alpha_samples[, dim(model$alpha_samples)[[2]]],
    num_obs = model$num_obs
  )

  carry_over <- c("metric", "type", "logz_list", "n_items", "n_particles",
                  "mcmc_kernel_app", "alpha_prop_sd", "lambda", "alpha_max",
                  "alpha", "aug_method", "leap_size")

  for(nm in carry_over) {
    eval(parse(text = paste0("ret$", nm, "<-model$", nm)))
  }
  ret$rankings <- rankings
  ret$num_obs <- nrow(rankings)

  class(ret) <- c("SMCMallowsNewUsers","SMCMallows")
  ret
}

