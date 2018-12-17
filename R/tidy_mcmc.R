tidy_mcmc <- function(fit){

  fit <- tidy_rho(fit)
  fit <- tidy_alpha(fit)
  fit <- tidy_cluster_assignment(fit)
  fit <- tidy_cluster_probabilities(fit)
  fit <- tidy_wcd(fit)
  fit <- tidy_augmented_data(fit)
  fit <- tidy_augmentation_acceptance(fit)
  fit <- tidy_error_probability(fit)

  return(fit)
}



tidy_rho <- function(fit){
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
                    levels = paste("Cluster", sort(unique(cluster))))

  iteration <- rep(seq(from = 1, to = rho_dims[[3]] * fit$rho_thinning, by = fit$rho_thinning),
                   each = rho_dims[[1]] * rho_dims[[2]])

  # Store the final rho as a tibble
  fit$rho <- dplyr::tibble(
    item = item,
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  return(fit)
}

tidy_alpha <- function(fit){
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
                    levels = paste("Cluster", sort(unique(cluster))))

  iteration <- rep(
    seq(from = 1, to = alpha_dims[[2]] * fit$alpha_jump, by = fit$alpha_jump),
    each = alpha_dims[[1]]
  )

  fit$alpha <- dplyr::tibble(
    cluster = cluster,
    iteration = iteration,
    value = value
  )

  return(fit)
}

tidy_cluster_assignment <- function(fit){

  # Tidy cluster assignment
  if(fit$save_clus){
    if(fit$n_clusters > 1){
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

    fit$cluster_assignment <- dplyr::tibble(
      assessor = assessor,
      iteration = iteration,
      value = value
    )
  } else {
    fit$cluster_assignment <- NULL
  }


  return(fit)
}

tidy_cluster_probabilities <- function(fit){
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
    seq(from = 1, to = clusprob_dims[[1]], by = 1),
    times = clusprob_dims[[2]]
  )

  cluster <- factor(paste("Cluster", cluster),
                    levels = paste("Cluster", sort(unique(cluster))))

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


tidy_wcd <- function(fit){
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
    cluster <- factor(paste("Cluster", cluster),
                      levels = paste("Cluster", sort(unique(cluster))))

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

tidy_augmented_data <- function(fit){
  # Tidy augmented data, or delete
  if(fit$save_aug){

    augdata_dims <- dim(fit$augmented_data)

    # Item1, Item2, ..., Item1, Item2, ..., Item1, Item2, ..., Item1, Item2
    # Assessor1, Assessor1, ..., Assessor2, Assessor2, ... Assessor1, Assessor1, ..., Assessor2, Assessor2
    # Iteration1, Iteration1, ..., Iteration1, Iteration1, ..., Iteration2, Iteration2, ... Iteration2, Iteration2
    value <- c(fit$augmented_data)

    item <- rep(fit$items, times = augdata_dims[[2]] * augdata_dims[[3]])
    item <- factor(item, levels = fit$items)

    assessor <- rep(seq(from = 1, to = augdata_dims[[2]], by = 1), each = augdata_dims[[1]],
                    times = augdata_dims[[3]])

    iteration <- rep(seq(from = 1, to = augdata_dims[[3]] * fit$aug_thinning, by = fit$aug_thinning),
                     each = augdata_dims[[1]] * augdata_dims[[2]])

    fit$augmented_data <- dplyr::tibble(
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


tidy_augmentation_acceptance <- function(fit){
  # Augmentation acceptance

  if(fit$any_missing || fit$augpair){
    fit$aug_acceptance <- dplyr::tibble(acceptance_rate = c(fit$aug_acceptance))
    fit$aug_acceptance <- dplyr::mutate(fit$aug_acceptance,
                                        assessor = dplyr::row_number())
    fit$aug_acceptance <- dplyr::select(fit$aug_acceptance,
                                        .data$assessor, .data$acceptance_rate)

  } else {
    fit$aug_acceptance <- NULL
  }
  return(fit)
}



tidy_error_probability <- function(fit){
  theta_length <- length(fit$theta)

  if(theta_length > 0){
    fit$theta <- dplyr::tibble(
      iteration = seq(from = 1, to = theta_length, by = 1),
      value = c(fit$theta)
    )
  } else {
    fit$theta <- NULL
  }



  return(fit)
}
