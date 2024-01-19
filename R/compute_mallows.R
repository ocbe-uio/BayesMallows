#' Preference Learning with the Mallows Rank Model
#'
#' @description Compute the posterior distributions of the parameters of the
#'   Bayesian Mallows Rank Model, given rankings or preferences stated by a set
#'   of assessors.
#'
#'   The `BayesMallows` package uses the following parametrization of the
#'   Mallows rank model \insertCite{mallows1957}{BayesMallows}:
#'
#'   \deqn{p(r|\alpha,\rho) = \frac{1}{Z_{n}(\alpha)} \exp\left\{\frac{-\alpha}{n}
#'   d(r,\rho)\right\}}
#'
#'   where \eqn{r} is a ranking, \eqn{\alpha} is a scale parameter, \eqn{\rho}
#'   is the latent consensus ranking, \eqn{Z_{n}(\alpha)} is the partition
#'   function (normalizing constant), and \eqn{d(r,\rho)} is a distance function
#'   measuring the distance between \eqn{r} and \eqn{\rho}. We refer to
#'   \insertCite{vitelli2018;textual}{BayesMallows} for further details of the Bayesian
#'   Mallows model.
#'
#'   `compute_mallows` always returns posterior distributions of the latent
#'   consensus ranking \eqn{\rho} and the scale parameter \eqn{\alpha}. Several
#'   distance measures are supported, and the preferences can take the form of
#'   complete or incomplete rankings, as well as pairwise preferences.
#'   `compute_mallows` can also compute mixtures of Mallows models, for
#'   clustering of assessors with similar preferences.
#'
#' @param data An object of class "BayesMallowsData" returned from
#'   [setup_rank_data()].
#'
#' @param model_options An object of class "BayesMallowsModelOptions" returned
#'   from [set_model_options()].
#'
#' @param compute_options An object of class "BayesMallowsComputeOptions"
#'   returned from [set_compute_options()].
#'
#' @param priors An object of class "BayesMallowsPriors" returned from
#'   [set_priors()].
#'
#' @param initial_values An object of class "BayesMallowsInitialValues" returned
#'   from [set_initial_values()].
#'
#' @param pfun_estimate Object returned from [estimate_partition_function()].
#'   Defaults to \code{NULL}, and will only be used for footrule, Spearman, or
#'   Ulam distances when the cardinalities are not available, cf.
#'   [get_cardinalities()].
#'
#' @param verbose Logical specifying whether to print out the progress of the
#'   Metropolis-Hastings algorithm. If `TRUE`, a notification is printed
#'   every 1000th iteration. Defaults to `FALSE`.
#'
#' @param cl Optional cluster returned from [parallel::makeCluster()]. If
#'   provided, chains will be run in parallel, one on each node of `cl`.
#'
#' @return An object of class BayesMallows.
#'
#' @references \insertAllCited{}
#'
#' @export
#' @importFrom rlang .data
#'
#' @family modeling
#'
#' @example /inst/examples/compute_mallows_example.R
#' @example /inst/examples/label_switching_example.R
#'
compute_mallows <- function(
    data,
    model_options = set_model_options(),
    compute_options = set_compute_options(),
    priors = set_priors(),
    initial_values = set_initial_values(),
    pfun_estimate = NULL,
    verbose = FALSE,
    cl = NULL) {
  validate_class(data, "BayesMallowsData")
  validate_class(model_options, "BayesMallowsModelOptions")
  validate_class(compute_options, "BayesMallowsComputeOptions")
  validate_class(priors, "BayesMallowsPriors")
  validate_class(initial_values, "BayesMallowsInitialValues")
  validate_preferences(data, model_options)
  validate_rankings(data)
  validate_initial_values(initial_values, data)

  pfun_values <- tryCatch(
    prepare_partition_function(model_options$metric, data$n_items),
    error = function(e) {
      if (is.null(pfun_estimate)) {
        stop(
          "Exact partition function not known. Please provide an ",
          "estimate in argument pfun_estimate."
        )
      } else {
        return(NULL)
      }
    }
  )

  if (is.null(cl)) {
    lapplyfun <- lapply
    chain_seq <- 1
  } else {
    lapplyfun <- prepare_cluster(cl, c(
      "data", "model_options", "compute_options", "priors", "initial_values",
      "pfun_values", "pfun_estimate", "verbose"
    ))
    chain_seq <- seq_along(cl)
  }

  fits <- lapplyfun(X = chain_seq, FUN = function(i) {
    if (length(initial_values$alpha_init) > 1) {
      initial_values$alpha_init <- initial_values$alpha_init[[i]]
    }
    run_mcmc(
      data = data,
      model_options = model_options,
      compute_options = compute_options,
      priors = priors,
      initial_values = initial_values,
      pfun_values = pfun_values,
      pfun_estimate = pfun_estimate,
      verbose = verbose
    )
  })

  fit <- tidy_mcmc(fits, data, model_options, compute_options)
  fit$pfun_values <- pfun_values
  fit$pfun_estimate <- pfun_estimate
  fit$priors <- priors
  class(fit) <- "BayesMallows"
  return(fit)
}
