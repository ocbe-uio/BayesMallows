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
  cluster_dims <- dim(fit$cluster_indicator)

  # Assessor1, Assessor2, ..., Assessor1, Assessor2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- c(fit$cluster_indicator)
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
  clusprob_dims <- dim(fit$cluster_probs)

  # Cluster1, Cluster2, ..., Cluster1, Cluster2
  # Iteration1, Iteration1, ..., Iteration2, Iteration2
  value <- c(fit$cluster_probs)

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

  return(fit)

}
