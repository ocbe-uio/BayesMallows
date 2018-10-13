# Then we do importance sampling for Spearman distance
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
system.time({estimates <- foreach(it = 15:22) %dopar% {
  BayesMallows::estimate_partition_function(alpha_vector = seq(from = 0, to = 10, by = .1),
                                            n_items = it, metric = "spearman", nmc = 1e3,
                                            degree = 9)
}})
stopCluster(cl)

save(estimates, file = "estimates.RData")
