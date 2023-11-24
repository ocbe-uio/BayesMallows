#' Set SMC compute options
#'
#' Sets the SMC compute options to be used in [update_mallows.BayesMallows()].
#'
#' @param n_particles Integer specifying the number of particles.
#' @param aug_method Augmentation proposal for use with missing data. One of
#'   "pseudo" and "uniform".
#' @param mcmc_steps Number of MCMC steps to be applied in the resample-move
#'   step.
#'
#' @return An object of class "SMCOptions".
#' @export
#'
#' @family preprocessing
set_smc_options <- function(n_particles = 1000, aug_method = "pseudo",
                            mcmc_steps = 5) {
  aug_method <- match.arg(aug_method, c("pseudo", "uniform"))

  validate_integer(n_particles)
  validate_integer(mcmc_steps)

  ret <- as.list(environment())
  class(ret) <- "SMCOptions"
  ret
}
