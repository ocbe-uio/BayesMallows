#' @title Set SMC compute options
#'
#' @description Sets the SMC compute options to be used in [update_mallows.BayesMallows()].
#'
#' @param n_particles Integer specifying the number of particles.
#' @param mcmc_steps Number of MCMC steps to be applied in the resample-move
#'   step.
#'
#' @return An object of class "SMCOptions".
#' @export
#'
#' @family preprocessing
set_smc_options <- function(
    n_particles = 1000,
    mcmc_steps = 5) {
  validate_integer(n_particles)
  validate_integer(mcmc_steps)

  ret <- as.list(environment())
  class(ret) <- "SMCOptions"
  ret
}
