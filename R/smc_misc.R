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
      from = burnin(model) + 1, to = ncol(model$alpha_samples),
      length.out = n_particles
    )
  )
  model$alpha_samples[1, thinned_inds, drop = TRUE]
}

extract_rho_init <- function(model, n_particles) {
  thinned_inds <- floor(
    seq(
      from = burnin(model) + 1, to = dim(model$rho_samples)[[3]],
      length.out = n_particles
    )
  )
  model$rho_samples[, 1, thinned_inds, drop = TRUE]
}

run_common_part <- function(
    data, new_data, model_options, smc_options, compute_options, priors,
    initial_values, pfun_list, model) {
  ret <- run_smc(
    data = data,
    new_data = list(new_data),
    model_options = model_options,
    smc_options = smc_options,
    compute_options = compute_options,
    priors = priors,
    initial_values = initial_values,
    pfun_values = pfun_list$pfun_values,
    pfun_estimate = pfun_list$pfun_estimate
  )

  ret <- c(ret, tidy_smc(ret, data$items))
  ret$model_options <- model_options
  ret$smc_options <- smc_options
  ret$compute_options <- compute_options
  class(ret$compute_options) <- "list"
  ret$priors <- priors
  ret$n_items <- model$n_items
  ret$n_clusters <- 1
  ret$data <- new_data
  ret$pfun_values <- pfun_list$pfun_values
  ret$pfun_estimate <- pfun_list$pfun_estimate
  ret$model_options$metric <- model_options$metric
  if(prod(dim(ret$augmented_rankings)) == 0) ret$augmented_rankings <- NULL
  ret$items <- data$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}

flush <- function(data) {
  data$rankings <- data$rankings[integer(), , drop = FALSE]
  data$n_assessors <- 0
  data$observation_frequency <- data$observation_frequency[integer()]
  data$consistent <- data$consistent[integer(), , drop = FALSE]
  data$user_ids <- data$user_ids[integer()]
  data$timepoint <- data$timepoint[integer()]
  data
}
