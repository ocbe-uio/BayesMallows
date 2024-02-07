#' @title Set SMC compute options
#'
#' @description Sets the SMC compute options to be used in
#'   [update_mallows.BayesMallows()].
#'
#' @param n_particles Integer specifying the number of particles.
#' @param mcmc_steps Number of MCMC steps to be applied in the resample-move
#'   step.
#' @param latent_sampling_lag Parameter specifying the number of timesteps to go
#'   back when resampling the latent ranks in the move step. See Section 6.2.3
#'   of \insertCite{Kantas2015}{BayesMallows} for details. The \eqn{L} in their
#'   notation corresponds to `latent_sampling_lag`. See more under Details.
#'   Defaults to `NA`, which means that all latent ranks from previous timesteps
#'   are resampled. If set to `0`, no move step is applied to the latent ranks.
#'
#' @return An object of class "SMCOptions".
#'
#' @section Lag parameter in move step:
#'
#'   The parameter `latent_sampling_lag` corresponds to \eqn{L} in
#'   \insertCite{Kantas2015}{BayesMallows}. Its use in this package is can be
#'   explained in terms of Algorithm 12 and Algorithm 13 in
#'   \insertCite{steinSequentialInferenceMallows2023}{BayesMallows}. In
#'   Algorithm 12, the relevant line is close to the bottom:
#'
#'   **for** \eqn{j = 1 : M_{t}} **do** \cr
#'   **M-H step:** update \eqn{\tilde{\mathbf{R}}_{j}^{(i)}} with proposal
#'   \eqn{\tilde{\mathbf{R}}_{j}' \sim q(\tilde{\mathbf{R}}_{j}^{(i)} |
#'   \mathbf{R}_{j}, \boldsymbol{\rho}_{t}^{(i)}, \alpha_{t}^{(i)})}.\cr
#'   **end**
#'
#'   Let \eqn{L} denote the value of `latent_sampling_lag`. With this parameter,
#'   we modify for algorithm so it becomes
#'
#'   **for** \eqn{j = M_{t-L+1} : M_{t}} **do** \cr
#'   **M-H step:** update \eqn{\tilde{\mathbf{R}}_{j}^{(i)}} with proposal
#'   \eqn{\tilde{\mathbf{R}}_{j}' \sim q(\tilde{\mathbf{R}}_{j}^{(i)} |
#'   \mathbf{R}_{j}, \boldsymbol{\rho}_{t}^{(i)}, \alpha_{t}^{(i)})}.\cr
#'   **end**
#'
#'   This means that with \eqn{L=0} no move step is performed on any latent
#'   ranks, whereas \eqn{L=1} means that the move step is only applied to the
#'   parameters entering at the given timestep. The default,
#'   `latent_sampling_lag = NA` means that \eqn{L=t} at each timestep, and
#'   hence all latent ranks are part of the move step at each timestep.
#'
#'
#' @export
#'
#' @family preprocessing
set_smc_options <- function(
    n_particles = 1000,
    mcmc_steps = 5,
    latent_sampling_lag = NA_integer_) {
  validate_integer(n_particles)
  validate_integer(mcmc_steps)

  ret <- as.list(environment())
  class(ret) <- "SMCOptions"
  ret
}
