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

# Then we add importance sampling estimates
# scp them with this command
# pscp oyss@abel.uio.no:/usit/abel/u1/oyss/BayesMallows/data-raw/importance_sampling/estimates.RData estimates.RData
load("./data-raw/importance_sampling/estimates.RData")

partition_function_data <- partition_function_data %>%
  bind_rows(estimates)

# Finally, keep the rows in pfd_old which are not in partition_function_data
partition_function_data <- pfd_old %>%
  anti_join(partition_function_data, by = c("n_items", "metric", "type")) %>%
  bind_rows(partition_function_data)

# Finally, save the fit as internal data
devtools::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
