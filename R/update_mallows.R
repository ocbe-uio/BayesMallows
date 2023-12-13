#' Update a Bayesian Mallows model with new users
#'
#' Update a Bayesian Mallows model estimated using the Metropolis-Hastings
#' algorithm in [compute_mallows()] using the sequential Monte Carlo algorithm
#' described in
#' \insertCite{steinSequentialInferenceMallows2023;textual}{BayesMallows}.
#'
#' @param model A model object of class "BayesMallows" returned from
#'   [compute_mallows()] or an object of class "SMCMallows" returned from
#'   this function.
#' @param new_data An object of class "BayesMallowsData" returned from
#'   [setup_rank_data()]. The object should contain the new data being provided.
#' @param model_options An object of class "BayesMallowsModelOptions" returned
#'   from [set_model_options()].
#' @param smc_options An object of class "SMCOptions" returned from
#'   [set_smc_options()].
#' @param compute_options An object of class "BayesMallowsComputeOptions"
#'   returned from [set_compute_options()].
#' @param priors An object of class "BayesMallowsPriors" returned from
#'   [set_priors()]. Defaults to the priors used in `model`.
#' @param ... Optional arguments. Currently not used.
#'
#' @return An updated model, of class "SMCMallows".
#' @export
#'
#' @family modeling
#'
#' @example /inst/examples/update_mallows_example.R
#'
update_mallows <- function(model, new_data, ...) {
  UseMethod("update_mallows")
}

#' @export
#' @rdname update_mallows
update_mallows.BayesMallows <- function(
    model, new_data,
    model_options = set_model_options(),
    smc_options = set_smc_options(),
    compute_options = set_compute_options(),
    priors = model$priors,
    ...) {
  if (is.null(model$burnin)) stop("Burnin must be set.")
  alpha_init <- extract_alpha_init(model, smc_options$n_particles)
  rho_init <- extract_rho_init(model, smc_options$n_particles)

  ret <- run_smc(
    data = new_data,
    new_data = new_data,
    model_options = model_options,
    smc_options = smc_options,
    compute_options = compute_options,
    priors = priors,
    initial_values = list(
      alpha_init = alpha_init, rho_init = rho_init,
      aug_init = NULL
    ),
    logz_list = model$logz_list
  )

  ret <- c(ret, tidy_smc(ret, model$items))

  ret$model_options <- model_options
  ret$smc_options <- smc_options
  ret$compute_options <- compute_options
  ret$priors <- priors
  ret$n_items <- model$n_items
  ret$burnin <- 0
  ret$n_clusters <- 1
  ret$data <- new_data
  ret$logz_list <- model$logz_list
  ret$metric <- model$metric
  ret$items <- model$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}


#' @export
#' @rdname update_mallows
update_mallows.SMCMallows <- function(model, new_data, ...) {
  if (!is.null(new_data$user_ids) && !is.null(model$data$user_ids)) {
    old_users <- setdiff(model$data$user_ids, new_data$user_ids)
    updated_users <- intersect(model$data$user_ids, new_data$user_ids)
    new_users <- setdiff(new_data$user_ids, model$data$user_ids)

    rankings <- rbind(
      model$data$rankings[old_users, , drop = FALSE],
      new_data$rankings[c(updated_users, new_users), , drop = FALSE]
    )

    user_ids <- c(old_users, updated_users, new_users)

    data <- setup_rank_data(rankings = rankings, user_ids = user_ids)
    new_data <- setup_rank_data(
      rankings = rankings[new_users, , drop = FALSE],
      user_ids = new_users
    )

    if (!is.null(model$augmented_rankings)) {
      consistent <- matrix(TRUE, nrow = nrow(rankings), ncol = model$smc_options$n_particles)
      for (uu in updated_users) {
        index <- which(rownames(rankings) == uu)
        to_compare <- as.numeric(stats::na.omit(rankings[index, ]))

        consistent[index, ] <- apply(model$augmented_rankings[, index, ], 2, function(ar) {
          all(ar[ar %in% to_compare] == to_compare)
        })
      }
      data$consistent <- consistent * 1L
    }
  } else {
    rankings <- rbind(model$data$rankings, new_data$rankings)
    data <- setup_rank_data(
      rankings = rankings,
      user_ids = seq_len(nrow(rankings))
    )
  }

  ret <- run_smc(
    data = data,
    new_data = new_data,
    model_options = model$model_options,
    smc_options = model$smc_options,
    compute_options = model$compute_options,
    priors = model$priors,
    initial_values = list(
      alpha_init = model$alpha_samples,
      rho_init = model$rho_samples,
      aug_init = model$augmented_rankings
    ),
    logz_list = model$logz_list
  )

  tidy_parameters <- tidy_smc(ret, model$items)
  model$alpha <- tidy_parameters$alpha
  model$rho <- tidy_parameters$rho
  model$augmented_rankings <- ret$augmented_rankings
  model$ESS <- ret$ESS
  model$data <- data

  class(model) <- c("SMCMallows", "BayesMallows")
  model
}


tidy_smc <- function(ret, items) {
  result <- list()
  result$alpha <- tidy_alpha(matrix(ret$alpha_samples, nrow = 1), 1, 1)

  rho_mat <- array(dim = c(dim(ret$rho_samples)[[1]], 1, dim(ret$rho_samples)[[2]]))
  rho_mat[, 1, ] <- ret$rho_samples
  result$rho <- tidy_rho(rho_mat, 1, 1, items)

  result
}

extract_alpha_init <- function(model, n_particles) {
  thinned_inds <- floor(
    seq(
      from = model$burnin + 1, to = ncol(model$alpha_samples),
      length.out = n_particles
    )
  )
  model$alpha_samples[1, thinned_inds, drop = TRUE]
}

extract_rho_init <- function(model, n_particles) {
  thinned_inds <- floor(
    seq(
      from = model$burnin + 1, to = dim(model$rho_samples)[[3]],
      length.out = n_particles
    )
  )
  model$rho_samples[, 1, thinned_inds, drop = TRUE]
}
