# We need dplyr and tidyr to do this nicely
library(dplyr)
library(tidyr)
library(purrr)

# First we load the cardinalities for the footrule
load("./data-raw/footrule_cardinalities.Rdata", verbose = TRUE)

# Then we create a tibble to hold the cardinalities
partition_function_data <- tibble(
  num_items = seq_along(seq2),
  metric = rep("footrule", length(seq2)),
  values = seq2
)
rm(seq2)

# Then we load the cardinalities for the Spearman distance
load("./data-raw/spearman_cardinalities.Rdata")

# We add these to the tibble
partition_function_data <- tibble(
  num_items = seq_along(seq2),
  metric = rep("spearman", length(seq2)),
  values = seq2
  ) %>%
  bind_rows(partition_function_data) %>%
  mutate(
    type = "cardinalities"
  )
rm(seq2)

# Next we add the importance sampling estimates to the same tibble

# Code for generating off-line importance sampling estimates for Spearman and footrule.

# tibble to hold the parameters
parameters <- bind_rows(
  crossing(num_items = 51:52, metric = "footrule", alpha = seq(from = 0.01, to = 10, length.out = 5), nmc = 1e2),
  crossing(num_items = 15:26, metric = "spearman", alpha = seq(from = 0.01, to = 20, length.out = 5), nmc = 1e2)
)

# Use pmap from purrr to estimate logZ for each combination of the parameters

importance_sampling_parameters <- parameters %>%
  pmap(function(num_items, metric, alpha, nmc) {
    tibble(
      num_items = num_items,
      metric = metric,
      alpha = alpha,
      logZ = as.numeric(BayesMallows:::compute_importance_sampling_estimate(alpha, num_items, metric, nmc))
      )
  }
  ) %>%
  map_dfr(.f = identity) %>%
  nest(-num_items, -metric) %>%
  mutate(
    values = map(data, ~ list(lm(.$logZ ~ poly(.$alpha, 4))$coefficients))
    ) %>%
  select(-data) %>%
  mutate(type = "importance_sampling")


# Finally, add these to the partition function data
partition_function_data <- partition_function_data %>%
  bind_rows(importance_sampling_parameters)


# Finally, save the fit as internal data
devtools::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
