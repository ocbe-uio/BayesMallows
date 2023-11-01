#' Preference Learning with the Mallows Rank Model
#'
#' @description Compute the posterior distributions of the parameters of the
#'   Bayesian Mallows Rank Model, given rankings or preferences stated by a set
#'   of assessors.
#'
#'   The \code{BayesMallows} package uses the following parametrization of the
#'   Mallows rank model \insertCite{mallows1957}{BayesMallows}:
#'
#'   \deqn{p(r|\alpha,\rho) = \frac{1}{Z_{n}(\alpha)} \exp\{\frac{-\alpha}{n}
#'   d(r,\rho)\}}
#'
#'   where \eqn{r} is a ranking, \eqn{\alpha} is a scale parameter, \eqn{\rho}
#'   is the latent consensus ranking, \eqn{Z_{n}(\alpha)} is the partition
#'   function (normalizing constant), and \eqn{d(r,\rho)} is a distance function
#'   measuring the distance between \eqn{r} and \eqn{\rho}. We refer to
#'   \insertCite{vitelli2018}{BayesMallows} for further details of the Bayesian
#'   Mallows model.
#'
#'   \code{compute_mallows} always returns posterior distributions of the latent
#'   consensus ranking \eqn{\rho} and the scale parameter \eqn{\alpha}. Several
#'   distance measures are supported, and the preferences can take the form of
#'   complete or incomplete rankings, as well as pairwise preferences.
#'   \code{compute_mallows} can also compute mixtures of Mallows models, for
#'   clustering of assessors with similar preferences.
#'
#' @param data An object of class \code{"BayesMallowsData"} returned from
#'   \code{\link{setup_rank_data}}.
#'
#' @param model An object of class \code{"BayesMallowsModelOptions"} returned
#'   from \code{\link{set_model_options}}.
#'
#' @param compute_options An object of class \code{"BayesMallowsComputeOptions"}
#'   returned from \code{\link{set_compute_options}}.
#'
#' @param priors An object of class \code{"BayesMallowsPriors"} returned from
#'   \code{\link{set_priors}}.
#'
#' @param init An object of class \code{"BayesMallowsInitialValues"} returned
#'   from \code{\link{set_initial_values}}.
#'
#' @param logz_estimate Estimate of the partition function, computed with
#'   \code{\link{estimate_partition_function}}. Defaults to \code{NULL}, which
#'   means means that \code{compute_mallows} attempts to compute the
#'   partition function by using either a pre-defined integer sequence or a
#'   pre-computed importance sampling estimate. When Cayley, Hamming, or
#'   Kendall distance is used, the partition function is efficiently
#'   computed analytically, and hence in these cases the argument is not needed,
#'   although it will still be used if provided.
#'
#' @param verbose Logical specifying whether to print out the progress of the
#'   Metropolis-Hastings algorithm. If \code{TRUE}, a notification is printed
#'   every 1000th iteration. Defaults to \code{FALSE}.
#'
#' @param seed Optional integer to be used as random number seed.
#'
#' @param cl Optional cluster returned from \code{parallel::makeCluster}. If
#'   provided, chains will be run in parallel, one on each node of \code{cl}.
#'
#'
#' @return A list of class BayesMallows.
#'
#' @seealso \code{\link{compute_mallows_mixtures}} for a function that computes
#'   separate Mallows models for varying numbers of clusters.
#'
#'
#'
#' @references \insertAllCited{}
#'
#' @export
#' @importFrom rlang .data
#'
#' @family modeling
#'
#' @example /inst/examples/compute_mallows_example.R
#'
compute_mallows <- function(
    data,
    model = set_model_options(),
    compute_options = set_compute_options(),
    priors = set_priors(),
    init = set_initial_values(),
    logz_estimate = NULL,
    verbose = FALSE,
    seed = NULL,
    cl = NULL) {

  validate_class(data, "BayesMallowsData")
  validate_class(model, "BayesMallowsModelOptions")
  validate_class(compute_options, "BayesMallowsComputeOptions")
  validate_class(priors, "BayesMallowsPriors")
  validate_class(init, "BayesMallowsInitialValues")
  validate_preferences(data, model)
  validate_initial_values(init, data)

  if (!is.null(seed)) set.seed(seed)

  data <- update_data(data, model)
  logz_list <- prepare_partition_function(logz_estimate, model$metric, data$n_items)

  if (compute_options$save_ind_clus) {
    abort <- readline(
      prompt = paste(
        compute_options$nmc, "csv files will be saved in your current working directory.",
        "Proceed? (yes/no): "
      )
    )
    if (tolower(abort) %in% c("n", "no")) stop()
  }

  if (is.null(cl)) {
    lapplyfun <- lapply
    chain_seq <- 1
  } else {
    parallel::clusterExport(
      cl = cl,
      varlist = c(
        "data", "model", "compute_options", "priors", "init",
        "logz_list", "verbose"
      ),
      envir = environment()
    )
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    lapplyfun <- function(X, FUN, ...) {
      parallel::parLapply(cl = cl, X = X, fun = FUN, ...)
    }
    chain_seq <- seq_along(cl)
  }
  # to extract one sample at a time. armadillo is column major, just like rankings
  fits <- lapplyfun(X = chain_seq, FUN = function(i) {
    if (length(init$alpha_init) > 1) {
      init$alpha_init <- init$alpha_init[[i]]
    }
    run_mcmc(
      data = data,
      model = model,
      compute_options = compute_options,
      priors = priors,
      init = init,
      logz_list = logz_list,
      verbose = verbose
    )
  })

  if (verbose) {
    print("Metropolis-Hastings algorithm completed. Post-processing data.")
  }

  fit <- tidy_mcmc(fits, data, model, compute_options)

  # Add class attribute
  class(fit) <- "BayesMallows"

  return(fit)
}

update_data <- function(data, model) {
  if (any(is.na(data$rankings))) {
    dn <- dimnames(data$rankings)
    data$rankings <- lapply(
      split(data$rankings, f = seq_len(nrow(data$rankings))),
      function(x) {
        if (sum(is.na(x)) == 1) x[is.na(x)] <- setdiff(seq_along(x), x)
        return(x)
      }
    )
    data$rankings <- do.call(rbind, data$rankings)
    dimnames(data$rankings) <- dn
  }

  data$constraints <- generate_constraints(data)
  if (is.null(data$obs_freq)) data$obs_freq <- rep(1, nrow(data$rankings))

  data
}
