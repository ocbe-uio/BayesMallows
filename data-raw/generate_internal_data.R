# This code generates the internal data sets needed for computing the partition functions.

# We need these packages to do this nicely
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

# Should be 1e6
nmc <- 1e5
# Should be at least 50
num_alphas <- 50

# list to hold the parameters
parameters <- crossing(
  num_items = 15:20,
  metric = "spearman",
  alpha = seq(from = 0.01, to = 20, length.out = num_alphas),
  nmc = nmc
  ) %>%
  split(seq(nrow(.)))

library(parallel)
cl <- makeCluster(7)
is_fit <- parLapply(cl = cl, X = parameters, fun = function(X){
  X[["logZ"]] <- BayesMallows::compute_importance_sampling_estimate(
    alpha_vector = X[["alpha"]],
    n = X[["num_items"]],
    metric = X[["metric"]],
    nmc = X[["nmc"]]
  )
  return(X)
})
stopCluster(cl)

# Now we turn the fits into a tibble again and estimate regressions per
# num_items and metric

is_fit <- bind_rows(is_fit) %>%
  map_dfr(.f = identity) %>%
  nest(-num_items, -metric) %>%
  mutate(
    values = map(data, ~ list(lm(.$logZ ~ poly(.$alpha, 4))$coefficients))
    ) %>%
  select(-data) %>%
  mutate(type = "importance_sampling")


# Finally, add these to the partition function data
partition_function_data <- partition_function_data %>%
  bind_rows(is_fit)


# Finally, save the fit as internal data
devtools::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
