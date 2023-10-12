tidy_mcmc <- function(fits, rho_thinning, rankings, alpha_jump,
                      n_clusters, nmc, aug_thinning, n_items) {
  fit <- list()

  # Add names of item
  if (!is.null(colnames(rankings))) {
    items <- colnames(rankings)
  } else {
    items <- paste("Item", seq(from = 1, to = n_items, by = 1))
  }

  fit$rho <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_rho(fits[[i]]$rho, i, rho_thinning, items)
  }))

  fit$alpha <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_alpha(fits[[i]]$alpha, i, alpha_jump)
  }))

  fit$cluster_assignment <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_cluster_assignment(fits[[i]]$cluster_assignment, i, n_clusters, fits[[i]]$n_assessors, nmc)
  }))

  fit$cluster_probs <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_cluster_probabilities(fits[[i]]$cluster_probs, i, n_clusters, nmc)
  }))

  fit$within_cluster_distance <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_wcd(fits[[i]]$within_cluster_distance, i)
  }))

  fit$augmented_data <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_augmented_data(fits[[i]]$augmented_data, i, items, aug_thinning)
  }))

  fit$aug_acceptance <- lapply(seq_along(fits), function(i) {
    tidy_augmentation_acceptance(
      fits[[i]]$aug_acceptance, i, fits[[i]]$any_missing, fits[[i]]$augpair
    )
  })

  fit$theta <- do.call(rbind, lapply(seq_along(fits), function(i) {
    tidy_error_probability(fits[[i]]$theta, i)
  }))

  fit$n_clusters <- n_clusters
  fit$items <- items
  fit$n_items <- n_items
  fit$n_assessors <- fits[[1]]$n_assessors
  fit$nmc <- nmc
  fit$alpha_acceptance <- rowMeans(matrix(
    vapply(fits, function(x) x$alpha_acceptance, numeric(n_clusters)),
    nrow = n_clusters
  ))
  fit$rho_acceptance <- rowMeans(matrix(
    vapply(fits, function(x) x$rho_acceptance, numeric(n_clusters)),
    nrow = n_clusters
  ))

  return(fit)
}



tidy_rho <- function(rho_mat, chain, rho_thinning, items) {
  # Tidy rho
  rho_dims <- dim(rho_mat)
  # Item1, Item2, Item3, ...., Item1, Item2, Item3
  # Cluster1, Cluster1, Cluster1, ..., Cluster2, Cluster2, Cluster2
  # Iteration1, Iteration1, ..., Iteration1, Iteration1, Iteration1, Iteration2
  value <- c(rho_mat)

  item <- rep(items, times = rho_dims[[2]] * rho_dims[[3]])
  item <- factor(item, levels = items)

  cluster <- rep(
    seq(from = 1, to = rho_dims[[2]], by = 1),
    each = rho_dims[[1]],
    times = rho_dims[[3]]
  )
  cluster <- factor(paste("Cluster", cluster),
    levels = paste("Cluster", sort(unique(cluster)))
  )

  iteration <- rep(seq(from = 1, to = rho_dims[[3]] * rho_thinning, by = rho_thinning),
    each = rho_dims[[1]] * rho_dims[[2]]
  )

  # Store the final rho as a dataframe
  data.frame(
    chain = factor(chain),
    item = item,
    cluster = cluster,
    iteration = iteration,
    value = value
  )
}

tidy_alpha <- function(alpha_mat, chain, alpha_jump) {
  # Tidy alpha
  alpha_dims <- dim(alpha_mat)
  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2


  cluster <- rep(
    seq(from = 1, to = alpha_dims[[1]], by = 1),
    times = alpha_dims[[2]]
  )

  cluster <- factor(paste("Cluster", cluster),
    levels = paste("Cluster", sort(unique(cluster)))
  )

  iteration <- rep(
    seq(from = 1, to = alpha_dims[[2]] * alpha_jump, by = alpha_jump),
    each = alpha_dims[[1]]
  )

  data.frame(
    chain = factor(chain),
    cluster = cluster,
    iteration = iteration,
    value = c(alpha_mat)
  )
}

tidy_cluster_assignment <- function(
    cluster_assignment, chain, n_clusters,
    n_assessors, nmc) {
  if (n_clusters > 1) {
    cluster_dims <- dim(cluster_assignment)
    value <- paste("Cluster", cluster_assignment)
  } else {
    cluster_dims <- c(n_assessors, nmc)
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

  data.frame(
    chain = factor(chain),
    assessor = assessor,
    iteration = iteration,
    value = value
  )
}

tidy_cluster_probabilities <- function(cluster_probs, chain, n_clusters, nmc) {
  # Tidy cluster probabilities
  if (n_clusters > 1) {
    clusprob_dims <- dim(cluster_probs)
    value <- c(cluster_probs)
  } else {
    clusprob_dims <- c(n_clusters, nmc)
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

  data.frame(
    chain = factor(chain),
    cluster = cluster,
    iteration = iteration,
    value = value
  )
}


tidy_wcd <- function(within_cluster_distance, chain) {
  # Tidy the within-cluster distances, or delete the empty matrix
  if (!is.null(within_cluster_distance)) {
    wcd_dims <- dim(within_cluster_distance)
    value <- c(within_cluster_distance)

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

    data.frame(
      chain = factor(chain),
      cluster = cluster,
      iteration = iteration,
      value = value
    )
  } else {
    NULL
  }
}

tidy_augmented_data <- function(augmented_data, chain, items, aug_thinning) {
  # Tidy augmented data, or delete
  if (!is.null(augmented_data) && prod(dim(augmented_data)) != 0) {
    augdata_dims <- dim(augmented_data)

    # Item1, Item2, ..., Item1, Item2, ..., Item1, Item2, ..., Item1, Item2
    # Assessor1, Assessor1, ..., Assessor2, Assessor2, ... Assessor1, Assessor1, ..., Assessor2, Assessor2
    # Iteration1, Iteration1, ..., Iteration1, Iteration1, ..., Iteration2, Iteration2, ... Iteration2, Iteration2
    value <- c(augmented_data)

    item <- rep(items, times = augdata_dims[[2]] * augdata_dims[[3]])
    item <- factor(item, levels = items)

    assessor <- rep(seq(from = 1, to = augdata_dims[[2]], by = 1),
      each = augdata_dims[[1]],
      times = augdata_dims[[3]]
    )

    iteration <- rep(seq(from = 1, to = augdata_dims[[3]] * aug_thinning, by = aug_thinning),
      each = augdata_dims[[1]] * augdata_dims[[2]]
    )

    data.frame(
      chain = factor(chain),
      iteration = iteration,
      assessor = assessor,
      item = item,
      value = value
    )
  } else {
    NULL
  }
}


tidy_augmentation_acceptance <- function(aug_acceptance, chain, any_missing, augpair) {
  # Augmentation acceptance

  if (any_missing || augpair) {
    aug_acceptance <- data.frame(acceptance_rate = c(aug_acceptance))
    aug_acceptance$assessor <- seq_len(nrow(aug_acceptance))
    aug_acceptance <- aug_acceptance[, c("assessor", "acceptance_rate")]
    aug_acceptance$chain <- chain
    aug_acceptance
  } else {
    NULL
  }
}

tidy_error_probability <- function(theta, chain) {
  theta_length <- length(theta)

  if (theta_length > 0) {
    data.frame(
      chain = factor(chain),
      iteration = seq(from = 1, to = theta_length, by = 1),
      value = c(theta)
    )
  } else {
    NULL
  }
}
