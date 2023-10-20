#' @title SMC-Mallows New Users
#' @description Function to perform resample-move SMC algorithm where we
#' receive new users with complete rankings at each time step. See Chapter 4
#' of \insertCite{steinSequentialInferenceMallows2023}{BayesMallows}
#'
#' @param R_obs Matrix containing the full set of observed rankings of size
#' n_assessors by n_items
#' @param type One of \code{"complete"}, \code{"partial"}, or
#' \code{"partial_alpha_fixed"}.
#' @param n_items Integer is the number of items in a ranking
#' @param metric A character string specifying the distance metric to use
#' in the Bayesian Mallows Model. Available options are \code{"footrule"},
#' \code{"spearman"}, \code{"cayley"}, \code{"hamming"}, \code{"kendall"}, and
#' \code{"ulam"}.
#' @param leap_size leap_size Integer specifying the step size of the
#' leap-and-shift proposal distribution
#' @param N Integer specifying the number of particles
#' @param Time Integer specifying the number of time steps in the SMC algorithm
#' @param logz_estimate Estimate of the partition function, computed with
#' \code{\link{estimate_partition_function}}.
#' @param cardinalities Cardinalities for exact evaluation of partition function,
#' returned from \code{\link{prepare_partition_function}}.
#' @param mcmc_kernel_app Integer value for the number of applications we
#' apply the MCMC move kernel
#' @param num_new_obs Integer value for the number of new observations
#' (complete rankings) for each time step
#' @param alpha_prop_sd Numeric value specifying the standard deviation of the
#'   lognormal proposal distribution used for \eqn{\alpha} in the
#'   Metropolis-Hastings algorithm. Defaults to \code{0.1}.
#' @param lambda Strictly positive numeric value specifying the rate parameter
#'   of the truncated exponential prior distribution of \eqn{\alpha}. Defaults
#'   to \code{0.1}. When \code{n_cluster > 1}, each mixture component
#'   \eqn{\alpha_{c}} has the same prior distribution.
#' @param alpha_max Maximum value of \code{alpha} in the truncated exponential
#'   prior distribution.
#' @param alpha A numeric value of the scale parameter which is known and fixed.
#' @param aug_method A character string specifying the approach for filling
#' in the missing data, options are "pseudolikelihood" or "random".
#' @param verbose Logical specifying whether to print out the progress of the
#' SMC-Mallows algorithm. Defaults to \code{FALSE}.
#'
#' @return a set of particles each containing a value of rho and alpha
#'
#' @export
#'
#' @example inst/examples/smc_mallows_new_users_complete_example.R
#'
#' @family modeling
#'
smc_mallows_new_users <- function(
  R_obs,
  type,
  n_items,
  N,
  Time,
  mcmc_kernel_app,
  num_new_obs,
  alpha_prop_sd = 0.5,
  lambda = 0.1,
  alpha_max = 1e6,
  alpha = 0,
  aug_method = "random",
  logz_estimate = NULL,
  cardinalities = NULL,
  verbose = FALSE,
  metric = "footrule",
  leap_size = 1
) {
  smc_mallows_new_users_cpp(
    R_obs,
    type,
    n_items,
    N,
    Time,
    mcmc_kernel_app,
    num_new_obs,
    alpha_prop_sd,
    lambda,
    alpha_max,
    alpha,
    aug_method,
    logz_estimate,
    cardinalities,
    verbose,
    metric,
    leap_size
  )
}
