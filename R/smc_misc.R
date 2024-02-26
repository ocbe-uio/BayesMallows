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

prepare_new_data <- function(model, new_data) {
  if (!is.null(new_data$user_ids) && !is.null(model$data$user_ids)) {
    old_users <- setdiff(model$data$user_ids, new_data$user_ids)
    updated_users <- intersect(model$data$user_ids, new_data$user_ids)
    new_users <- setdiff(new_data$user_ids, model$data$user_ids)

    rankings <- rbind(
      model$data$rankings[old_users, , drop = FALSE],
      new_data$rankings[updated_users, , drop = FALSE]
    )
    data <- setup_rank_data(
      rankings = rankings, user_ids = c(old_users, updated_users)
    )
    new_data <- setup_rank_data(
      rankings = new_data$rankings[new_users, , drop = FALSE],
      user_ids = new_users
    )

    if (!is.null(model$augmented_rankings)) {
      consistent <- matrix(
        TRUE,
        nrow = nrow(rankings), ncol = model$smc_options$n_particles
      )
      for (uu in updated_users) {
        index <- which(rownames(rankings) == uu)
        to_compare <- as.numeric(stats::na.omit(rankings[index, ]))

        consistent[index, ] <- apply(model$augmented_rankings[, index, ], 2, function(ar) {
          all(ar[ar %in% to_compare] == to_compare)
        })
      }
      data$consistent <- consistent * 1L
      new_data$consistent <- matrix(1L,
        nrow = nrow(new_data$rankings),
        ncol = ncol(data$consistent)
      )
    }
  } else {
    data <- model$data
  }
  list(data = data, new_data = new_data)
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
  ret$priors <- priors
  ret$n_items <- model$n_items
  ret$n_clusters <- 1
  ret$data <- new_data
  ret$pfun_values <- pfun_list$pfun_values
  ret$pfun_estimate <- pfun_list$pfun_estimate
  ret$model_options$metric <- model_options$metric
  ret$items <- data$items
  class(ret) <- c("SMCMallows", "BayesMallows")
  ret
}

flush <- function(data) {
  if(!is.null(data$rankings)) data$rankings <- data$rankings[integer(), , drop = FALSE]
  data$n_assessors <- 0
  data$consistent <- data$consistent[integer(), , drop = FALSE]
  if(!is.null(data$user_ids)) data$user_ids <- data$user_ids[integer()]
  data$timepoint <- data$timepoint[integer()]
  data
}


join_data <- function(data, new_data) {
  data$n_assessors <- data$n_assessors + new_data$n_assessors
  data$consistent <- rbind(data$consistent, new_data$consistent)
  data$rankings <- rbind(data$rankings, new_data$rankings)
  if(!is.null(data$user_ids)) data$user_ids <- c(data$user_ids, new_data$user_ids)
  data$observation_frequency <- c(data$observation_frequency, new_data$observation_frequency)
  data$timepoint <- c(data$timepoint, new_data$timepoint)
  data
}
