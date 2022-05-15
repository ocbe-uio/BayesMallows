#' @describeIn smc_mallows_new_users Deprecated function for
#' \code{smc_mallows_new_users_partial}.
#' @export
smc_mallows_new_users_partial <- function(
    R_obs, n_items, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app,
    num_new_obs, alpha_prop_sd, lambda, alpha_max, aug_method, verbose = FALSE
) {
  .Deprecated("smc_mallows_new_users")
  smc_mallows_new_users(
    R_obs, "partial", n_items, metric, leap_size, N, Time, mcmc_kernel_app,
    num_new_obs, alpha_prop_sd, lambda, alpha_max, 0, aug_method, logz_estimate,
    verbose
  )
}

#' @describeIn smc_mallows_new_users Deprecated function for
#' \code{smc_mallows_new_users_complete}.
#' @export
smc_mallows_new_users_complete <- function(
  R_obs, n_items, metric, leap_size, N, Time, mcmc_kernel_app, num_new_obs,
  alpha_prop_sd, lambda, alpha_max, logz_estimate, verbose = FALSE
) {
  .Deprecated("smc_mallows_new_users")
  smc_mallows_new_users(
    R_obs, "complete", n_items, metric, leap_size, N, Time, mcmc_kernel_app,
    num_new_obs, alpha_prop_sd, lambda, alpha_max, 0, "random",
    logz_estimate, verbose
  )
}

#' @describeIn smc_mallows_new_users Deprecated function for
#' \code{smc_mallows_new_users_partial_alpha_fixed}.
#' @export
smc_mallows_new_users_partial_alpha_fixed <- function(
    R_obs, n_items, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app,
    num_new_obs, aug_method, alpha) {
  .Deprecated("smc_mallows_new_users")
  smc_mallows_new_users(
    R_obs, "partial_alpha_fixed", n_items, metric, leap_size, N, Time,
    mcmc_kernel_app, num_new_obs, 1, 1, 1, alpha, aug_method, logz_estimate,
    FALSE)
}

#' @describeIn smc_mallows_new_item_rank Deprecated function for
#' \code{smc_mallows_new_item_rank_alpha_fixed}.
#' @export
smc_mallows_new_item_rank_alpha_fixed <- function(
  alpha, n_items, R_obs, metric, leap_size, N, Time, logz_estimate,
  mcmc_kernel_app, alpha_prop_sd, lambda, alpha_max, aug_method,
  verbose = FALSE
){
  .Deprecated("smc_mallows_new_item_rank")
  smc_mallows_new_item_rank(
    n_items, R_obs, metric, leap_size, N, Time, logz_estimate, mcmc_kernel_app,
    alpha, alpha_prop_sd, lambda, alpha_max, aug_method, verbose, TRUE)
}
