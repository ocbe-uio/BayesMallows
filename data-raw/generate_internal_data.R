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



# Finally, save the fit as internal data
devtools::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
