# Then we do importance sampling for Spearman distance
n_items <- seq(from = 30L, to = 40L, by = 1L)
num_workers <- as.integer(Sys.getenv("SLURM_NTASKS")) - 1
cl <- parallel::makeCluster(num_workers)

doParallel::registerDoParallel(cl)
system.time({
  estimates <- foreach::foreach(it = n_items) %dopar% {
    BayesMallows::estimate_partition_function(
      alpha_vector = seq(from = 0, to = 20, by = .1),
      n_items = it, metric = "spearman", nmc = 1e6,
      degree = 9
    )
  }
})

parallel::stopCluster(cl)

estimates <- dplyr::tibble(
  n_items = n_items,
  values = estimates
) %>%
  dplyr::mutate(
    metric = "spearman",
    type = "importance_sampling",
    message = "Partition function estimated with importance sampling for alpha between 0 and 20."
  )

save(estimates, file = "estimates.RData")
