# Then we do importance sampling for Spearman distance
library(parallel)
n_items <- seq(from = 15L, to = 22L, by = 1L)
num_workers <- as.integer(Sys.getenv("SLURM_NTASKS")) - 1
myCluster <- makeCluster(num_workers, type = "MPI")

estimates <- parLapply(myCluster, n_items, function(n_items) {
  BayesMallows::estimate_partition_function(method = "importance_sampling",
                                            alpha_vector = seq(from = 0, to = 20, by = 0.1),
                                            n_items = n_items, metric = "spearman", nmc = 1e3,
                                            degree = 9)
})

stopCluster(myCluster)
Rmpi::mpi.finalize()

library(dplyr)
estimates <- tibble(
  n_items = n_items,
  values = estimates
) %>%
  mutate(
    metric = "spearman",
    type = "importance_sampling",
    message = "Partition function estimated with importance sampling for alpha between 0 and 20."
  )

save(estimates, file = "estimates.RData")
