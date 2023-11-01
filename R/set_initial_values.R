#' Set initial values of scale parameter and modal ranking
#'
#' @param rho_init Numeric vector specifying the initial value of the latent
#'   consensus ranking \eqn{\rho}. Defaults to NULL, which means that the
#'   initial value is set randomly. If `rho_init` is provided when
#'   `n_clusters > 1`, each mixture component \eqn{\rho_{c}} gets the same
#'   initial value.
#'
#'
#' @param alpha_init Numeric value specifying the initial value of the scale
#'   parameter \eqn{\alpha}. Defaults to `1`. When `n_clusters > 1`,
#'   each mixture component \eqn{\alpha_{c}} gets the same initial value. When
#'   chains are run in parallel, by providing an argument `cl = cl`, then
#'   `alpha_init` can be a vector of of length `length(cl)`, each
#'   element of which becomes an initial value for the given chain.
#'
#'
#' @return An object of class `"BayesMallowsInitialValues"`, to be
#'   provided to the `initial_values` argument of [compute_mallows()].
#' @export
#'
#' @family options
#'
set_initial_values <- function(rho_init = NULL, alpha_init = 1) {

  if (!is.null(rho_init)) {
    if (!validate_permutation(rho_init)) stop("rho_init must be a proper permutation")
    if (!(sum(is.na(rho_init)) == 0)) stop("rho_init cannot have missing values")
    rho_init <- matrix(rho_init, ncol = 1)
  }

  ret <- as.list(environment())
  class(ret) <- "BayesMallowsInitialValues"
  ret

}
