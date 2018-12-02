# This code generates the internal data sets needed for computing the partition functions.
# Start with the current
pfd_old <- BayesMallows:::partition_function_data

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
    type = "cardinalities",
    message = "Using exact partition function."
  )
rm(seq2)

# Then we compute the integer sequences for Ulam distance
library(doParallel)
cl <- makeCluster(8)
registerDoParallel(cl)
n_items <- 50

seq_ulam <- foreach(ni = 1:n_items) %dopar% {
  Z <- numeric(ni)
  for(d in 0:(ni - 1)){
    Z[d + 1] <- PerMallows::count.perms(perm.length = ni, dist.value = d, dist.name = "ulam")
  }
  Z
}


# We add these to the tibble
partition_function_data <- tibble(
  n_items = 1:n_items,
  metric = rep("ulam", length(seq_ulam)),
  values = seq_ulam
) %>%
  bind_rows(partition_function_data) %>%
  mutate(
    type = "cardinalities",
    message = "Using exact partition function."
  )
rm(seq2)


# Then we add importance sampling estimates
load("./data-raw/importance_sampling/estimates.RData")

partition_function_data <- partition_function_data %>%
  bind_rows(estimates)




# Finally, keep the rows in pfd_old which are not in partition_function_data
partition_function_data <- pfd_old %>%
  anti_join(partition_function_data, by = c("n_items", "metric", "type")) %>%
  bind_rows(partition_function_data)

# Finally, save the fit as internal data
usethis::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
