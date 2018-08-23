tidy_mcmc <- function(fit){
  # Tidy rho
  rho_dims <- dim(fit$rho)
  # Item1, Item2, Item3, ...., Item1, Item2, Item3
  # Cluster1, Cluster1, Cluster1, ..., Cluster2, Cluster2, Cluster2
  # Iteration1, Iteration1, ..., Iteration1, Iteration1, Iteration1, Iteration2
  value <- c(fit$rho)

  item <- rep(dimnames(fit$rho)[[1]], times = rho_dims[[2]] * rho_dims[[3]])

  cluster <- rep(
    paste("Cluster", seq(from = 1, to = rho_dims[[2]], by = 1)),
    each = rho_dims[[1]],
    times = rho_dims[[3]]
    )

  iteration <- rep(seq(from = 1, to = rho_dims[[3]], by = 1),
                   each = rho_dims[[1]] * rho_dims[[2]])

  # Store the final rho as a tibble
  fit$rho <- dplyr::tibble(
    item = item,
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  # Tidy alpha
  alpha_dims <- dim(fit$alpha)
  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- c(fit$alpha)

  cluster <- rep(
    paste("Cluster", seq(from = 1, to = alpha_dims[[1]], by = 1)),
    times = alpha_dims[[2]]
    )

  iteration <- rep(
    seq(from = 1, to = alpha_dims[[2]], by = 1),
    each = alpha_dims[[1]]
    )

  fit$alpha <- dplyr::tibble(
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  # Tidy cluster indicator
  if(fit$n_clusters > 1){
    cluster_dims <- dim(fit$cluster_indicator)
    value <- paste("Cluster", c(fit$cluster_indicator))
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

  fit$cluster_indicator <- dplyr::tibble(
    assessor = assessor,
    iteration = iteration,
    value = value
  )

  # Tidy cluster probabilities
  if(fit$n_clusters > 1){
    clusprob_dims <- dim(fit$cluster_probs)
    value <- c(fit$cluster_probs)
  } else {
    clusprob_dims <- c(fit$n_clusters, fit$nmc)
    value <- rep(1, times = prod(clusprob_dims))
  }


  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  cluster <- rep(
    paste("Cluster", seq(from = 1, to = clusprob_dims[[1]], by = 1)),
    times = clusprob_dims[[2]]
  )

  iteration <- rep(
    seq(from = 1, to = clusprob_dims[[2]], by = 1),
    each = clusprob_dims[[1]]
  )

  fit$cluster_probs <- dplyr::tibble(
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  # Tidy the within-cluster distances, or delete the empty matrix
  if(fit$include_wcd){
    wcd_dims <- dim(fit$within_cluster_distance)
    value <- c(fit$within_cluster_distance)

    # Cluster1, Cluster2, ..., Cluster1, Cluster2
    # Iteration1, Iteration1, ..., Iteration2, Iteration2

    cluster <- rep(
      paste("Cluster", seq(from = 1, to = wcd_dims[[1]], by = 1)),
      times = wcd_dims[[2]]
    )

    iteration <- rep(
      seq(from = 1, to = wcd_dims[[2]], by = 1),
      each = wcd_dims[[1]]
    )

    fit$within_cluster_distance <- dplyr::tibble(
      cluster = cluster,
      iteration = iteration,
      value = value
    )


  } else {
    fit$within_cluster_distance <- NULL
  }


  return(fit)

}
