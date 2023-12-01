#' Specify options for Metropolis-Hastings algorithms
#'
#' @param nmc Integer specifying the number of iteration of the
#'   Metropolis-Hastings algorithm to run. Defaults to `2000`. See
#'   [assess_convergence()] for tools to check convergence of the Markov chain.
#'
#' @param burnin Integer defining the number of samples to discard. Defaults to
#'   `NULL`, which means that burn-in is not set.
#'
#' @param alpha_prop_sd Numeric value specifying the \eqn{\sigma}
#'   parameter of the lognormal proposal distribution used for \eqn{\alpha} in
#'   the Metropolis-Hastings algorithm. The logarithm of the proposed samples
#'   will have standard deviation given by `alpha_prop_sd`. Defaults to `0.1`.
#'
#' @param leap_size Integer specifying the step size of the leap-and-shift
#'   proposal distribution used for proposed new latent ranks \eqn{rho}.
#'   Defaults to 1.
#'
#' @param swap_leap Integer specifying the step size of the swap proposal
#'   used for augmentation with pairwise preference data. Defaults to 1.
#'
#' @param alpha_jump Integer specifying how many times to sample \eqn{\rho}
#'   between each sampling of \eqn{\alpha}. In other words, how many times to
#'   jump over \eqn{\alpha} while sampling \eqn{\rho}, and possibly other
#'   parameters like augmented ranks \eqn{\tilde{R}} or cluster assignments
#'   \eqn{z}. Setting `alpha_jump` to a high number can speed up computation
#'   time, by reducing the number of times the partition function for the
#'   Mallows model needs to be computed. Defaults to `1`.
#'
#' @param aug_thinning Integer specifying the thinning for saving augmented
#'   data. Only used when `save_aug = TRUE`. Defaults to `1`.
#'
#' @param clus_thinning Integer specifying the thinning to be applied to cluster
#'   assignments and cluster probabilities. Defaults to `1`.
#'
#' @param rho_thinning Integer specifying the thinning of `rho` to be performed
#'   in the Metropolis- Hastings algorithm. Defaults to `1`. `compute_mallows`
#'   save every `rho_thinning`th value of \eqn{\rho}.
#'
#' @param include_wcd Logical indicating whether to store the within-cluster
#'   distances computed during the Metropolis-Hastings algorithm. Defaults to
#'   `FALSE`. Setting `include_wcd = TRUE` is useful when deciding the number of
#'   mixture components to include, and is required by [plot_elbow()].
#'
#' @param save_aug Logical specifying whether or not to save the augmented
#'   rankings every `aug_thinning`th iteration, for the case of missing data or
#'   pairwise preferences. Defaults to `FALSE`. Saving augmented data is useful
#'   for predicting the rankings each assessor would give to the items not yet
#'   ranked, and is required by [plot_top_k()].
#'
#' @param save_ind_clus Whether or not to save the individual cluster
#'   probabilities in each step. This results in csv files `cluster_probs1.csv`,
#'   `cluster_probs2.csv`, ..., being saved in the calling directory. This
#'   option may slow down the code considerably, but is necessary for detecting
#'   label switching using Stephen's algorithm.
#'
#'
#'
#' @return An object of class `"BayesMallowsComputeOptions"`, to be provided in
#'   the `compute_options` argument to [compute_mallows()].
#' @export
#'
#' @family preprocessing
#'
set_compute_options <- function(
    nmc = 2000,
    burnin = NULL,
    alpha_prop_sd = 0.1,
    leap_size = 1,
    swap_leap = 1,
    alpha_jump = 1,
    aug_thinning = 1,
    clus_thinning = 1,
    rho_thinning = 1,
    include_wcd = FALSE,
    save_aug = FALSE,
    save_ind_clus = FALSE) {
  validate_integer(nmc)
  validate_positive(alpha_prop_sd)
  validate_integer(leap_size)
  validate_integer(swap_leap)
  validate_integer(alpha_jump)
  validate_integer(aug_thinning)
  validate_integer(clus_thinning)
  validate_integer(rho_thinning)
  validate_logical(include_wcd)
  validate_logical(save_aug)
  validate_logical(save_ind_clus)

  check_larger(nmc, alpha_jump)
  check_larger(nmc, aug_thinning)
  check_larger(nmc, clus_thinning)
  check_larger(nmc, rho_thinning)

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsComputeOptions"
  ret
}


prompt_save_files <- function(compute_options) {
  if (compute_options$save_ind_clus) {
    abort <- readline(
      prompt = paste(
        compute_options$nmc, "csv files will be saved in your current working directory.",
        "Proceed? (yes/no): "
      )
    )
    if (tolower(abort) %in% c("n", "no")) stop()
  }
}
