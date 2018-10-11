# This code generates the internal data sets needed for computing the partition functions.

# We need these packages to do this nicely
library(dplyr)

# First we load the cardinalities for the footrule
load("./data-raw/footrule_cardinalities.Rdata", verbose = TRUE)

# Then we create a tibble to hold the cardinalities
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("footrule", length(seq2)),
  values = seq2
)
rm(seq2)

# Then we load the cardinalities for the Spearman distance
load("./data-raw/spearman_cardinalities.Rdata")

# We add these to the tibble
partition_function_data <- tibble(
  n_items = seq_along(seq2),
  metric = rep("spearman", length(seq2)),
  values = seq2
  ) %>%
  bind_rows(partition_function_data) %>%
  mutate(
    type = "cardinalities"
  )
rm(seq2)

# Then we do importance sampling for Spearman distance
library(doParallel)
cl <- makeCluster(6)
registerDoParallel(cl)
system.time({estimates <- foreach(it = 15:22) %dopar% {
  BayesMallows::estimate_partition_function(alpha_vector = seq(from = 0, to = 10, by = .1),
                                            n_items = it, metric = "spearman", nmc = 1e6,
                                            degree = 9)
}})
stopCluster(cl)

# Finally, save the fit as internal data
devtools::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
