library(tidyverse)
library(BayesMallows)

set.seed(1)
cluster1 <- sample_mallows(rho0 = 1:5, alpha0 = 5, n_samples = 20)
cluster2 <- sample_mallows(rho0 = c(1, 3, 4, 2, 5), alpha0 = 4, n_samples = 20)
cluster3 <- sample_mallows(rho0 = c(1, 2, 5, 4, 3), alpha0 = 3, n_samples = 20)

cluster_data <- rbind(cluster1, cluster2, cluster3)

usethis::use_data(cluster_data, overwrite = TRUE)
