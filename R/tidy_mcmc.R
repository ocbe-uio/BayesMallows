tidy_mcmc <- function(fit, chain = 1) {
  fit <- tidy_rho(fit, chain)
  fit <- tidy_alpha(fit, chain)
  fit <- tidy_cluster_assignment(fit, chain)
  fit <- tidy_cluster_probabilities(fit, chain)
  fit <- tidy_wcd(fit, chain)
  fit <- tidy_augmented_data(fit, chain)
  fit <- tidy_augmentation_acceptance(fit, chain)
  fit <- tidy_error_probability(fit, chain)

  return(fit)
}



tidy_rho <- function(fit, chain) {
  # Tidy rho
  rho_dims <- dim(fit$rho)
  # Item1, Item2, Item3, ...., Item1, Item2, Item3
  # Cluster1, Cluster1, Cluster1, ..., Cluster2, Cluster2, Cluster2
  # Iteration1, Iteration1, ..., Iteration1, Iteration1, Iteration1, Iteration2
  value <- c(fit$rho)

  item <- rep(fit$items, times = rho_dims[[2]] * rho_dims[[3]])
  item <- factor(item, levels = fit$items)

  cluster <- rep(
    seq(from = 1, to = rho_dims[[2]], by = 1),
    each = rho_dims[[1]],
    times = rho_dims[[3]]
  )
  cluster <- factor(paste("Cluster", cluster),
    levels = paste("Cluster", sort(unique(cluster)))
  )

  iteration <- rep(seq(from = 1, to = rho_dims[[3]] * fit$rho_thinning, by = fit$rho_thinning),
    each = rho_dims[[1]] * rho_dims[[2]]
  )

  # Store the final rho as a dataframe
  fit$rho <- data.frame(
    chain = factor(chain),
    item = item,
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  return(fit)
}

tidy_alpha <- function(fit, chain) {
  # Tidy alpha
  alpha_dims <- dim(fit$alpha)
  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- c(fit$alpha)

  cluster <- rep(
    seq(from = 1, to = alpha_dims[[1]], by = 1),
    times = alpha_dims[[2]]
  )

  cluster <- factor(paste("Cluster", cluster),
    levels = paste("Cluster", sort(unique(cluster)))
  )

  iteration <- rep(
    seq(from = 1, to = alpha_dims[[2]] * fit$alpha_jump, by = fit$alpha_jump),
    each = alpha_dims[[1]]
  )

  fit$alpha <- data.frame(
    chain = factor(chain),
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  return(fit)
}

tidy_cluster_assignment <- function(fit, chain) {
  # Tidy cluster assignment
  if (fit$save_clus) {
    if (fit$n_clusters > 1) {
      cluster_dims <- dim(fit$cluster_assignment)
      value <- paste("Cluster", c(fit$cluster_assignment))
    } else {
      cluster_dims <- c(fit$n_assessors, fit$nmc)
      value <- paste("Cluster", rep(1, prod(cluster_dims)))
    }


    # Assessor1, Assessor2, ..., Assessor1, Assessor2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2

    assessor <- rep(
      seq(from = 1, to = cluster_dims[[1]], by = 1),
      times = cluster_dims[[2]]
    )
    iteration <- rep(
      seq(from = 1, to = cluster_dims[[2]], by = 1),
      each = cluster_dims[[1]]
    )

    fit$cluster_assignment <- data.frame(
      chain = factor(chain),
      assessor = assessor,
      iteration = iteration,
      value = value
    )
  } else {
    fit$cluster_assignment <- NULL
  }


  return(fit)
}

tidy_cluster_probabilities <- function(fit, chain) {
  # Tidy cluster probabilities
  if (fit$n_clusters > 1) {
    clusprob_dims <- dim(fit$cluster_probs)
    value <- c(fit$cluster_probs)
  } else {
    clusprob_dims <- c(fit$n_clusters, fit$nmc)
    value <- rep(1, times = prod(clusprob_dims))
  }


  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  cluster <- rep(
    seq(from = 1, to = clusprob_dims[[1]], by = 1),
    times = clusprob_dims[[2]]
  )

  cluster <- factor(paste("Cluster", cluster),
    levels = paste("Cluster", sort(unique(cluster)))
  )

  iteration <- rep(
    seq(from = 1, to = clusprob_dims[[2]], by = 1),
    each = clusprob_dims[[1]]
  )

  fit$cluster_probs <- data.frame(
    chain = factor(chain),
    cluster = cluster,
    iteration = iteration,
    value = value
  )
  return(fit)
}


tidy_wcd <- function(fit, chain) {
  # Tidy the within-cluster distances, or delete the empty matrix
  if (fit$include_wcd) {
    wcd_dims <- dim(fit$within_cluster_distance)
    value <- c(fit$within_cluster_distance)

    # Cluster1, Cluster2, ..., Cluster1, Cluster2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2

    cluster <- rep(
      paste("Cluster", seq(from = 1, to = wcd_dims[[1]], by = 1)),
      times = wcd_dims[[2]]
    )
    cluster <- factor(paste("Cluster", cluster),
      levels = paste("Cluster", sort(unique(cluster)))
    )

    iteration <- rep(
      seq(from = 1, to = wcd_dims[[2]], by = 1),
      each = wcd_dims[[1]]
    )

    fit$within_cluster_distance <- data.frame(
      chain = factor(chain),
      cluster = cluster,
      iteration = iteration,
      value = value
    )
  } else {
    fit$within_cluster_distance <- NULL
  }
  return(fit)
}

tidy_augmented_data <- function(fit, chain) {
  # Tidy augmented data, or delete
  if (fit$save_aug) {
    augdata_dims <- dim(fit$augmented_data)

    # Item1, Item2, ..., Item1, Item2, ..., Item1, Item2, ..., Item1, Item2
    # Assessor1, Assessor1, ..., Assessor2, Assessor2, ... Assessor1, Assessor1, ..., Assessor2, Assessor2
    # Iteration1, Iteration1, ..., Iteration1, Iteration1, ..., Iteration2, Iteration2, ... Iteration2, Iteration2
    value <- c(fit$augmented_data)

    item <- rep(fit$items, times = augdata_dims[[2]] * augdata_dims[[3]])
    item <- factor(item, levels = fit$items)

    assessor <- rep(seq(from = 1, to = augdata_dims[[2]], by = 1),
      each = augdata_dims[[1]],
      times = augdata_dims[[3]]
    )

    iteration <- rep(seq(from = 1, to = augdata_dims[[3]] * fit$aug_thinning, by = fit$aug_thinning),
      each = augdata_dims[[1]] * augdata_dims[[2]]
    )

    fit$augmented_data <- data.frame(
      chain = factor(chain),
      iteration = iteration,
      assessor = assessor,
      item = item,
      value = value
    )
  } else {
    fit$augmented_data <- NULL
  }

  return(fit)
}


tidy_augmentation_acceptance <- function(fit, chain) {
  # Augmentation acceptance

  if (fit$any_missing || fit$augpair) {
    fit$aug_acceptance <- data.frame(acceptance_rate = c(fit$aug_acceptance))
    fit$aug_acceptance$assessor <- seq_len(nrow(fit$aug_acceptance))
    fit$aug_acceptance <- fit$aug_acceptance[, c("assessor", "acceptance_rate")]
    fit$aug_acceptance$chain <- chain
  } else {
    fit$aug_acceptance <- NULL
  }
  return(fit)
}



tidy_error_probability <- function(fit, chain) {
  theta_length <- length(fit$theta)

  if (theta_length > 0) {
    fit$theta <- data.frame(
      chain = factor(chain),
      iteration = seq(from = 1, to = theta_length, by = 1),
      value = c(fit$theta)
    )
  } else {
    fit$theta <- NULL
  }



  return(fit)
}
