# This code generates the internal data sets needed for computing the partition functions.
# Start with the current
pfd_old <- partition_function_data
# First we load the cardinalities for the footrule
load("./data-raw/footrule_cardinalities.Rdata", verbose = TRUE)

# Then we create a tibble to hold the cardinalities for footrule
cardinalities_footrule_df <- data.frame(
  n_items = seq_along(seq2),
  metric = rep("footrule", length(seq2))
)
cardinalities_footrule_df$values <- seq2
rm(seq2)

# Then we load the cardinalities for the Spearman distance
load("./data-raw/spearman_cardinalities.Rdata")

# We add these to the tibble
cardinalities_spearman_df <- data.frame(
  n_items = seq_along(seq2),
  metric = rep("spearman", length(seq2))
)
cardinalities_spearman_df$values <- seq2
rm(seq2)

# Then we compute the integer sequences for Ulam distance
n_items <- 30

seq_ulam <- list(n_items)
for (ni in 1:n_items) {
  seq_ulam[[ni]] <- numeric(ni)
  for (d in 0:(ni - 1)) {
    seq_ulam[[ni]][d + 1] <- PerMallows::count.perms(perm.length = ni, dist.value = d, dist.name = "ulam")
  }
}

cardinalities_ulam_df <- data.frame(
  n_items = 1:n_items,
  metric = rep("ulam", length(seq_ulam))
)
cardinalities_ulam_df$values <- seq_ulam
rm(seq_ulam)

cardinalities_df <- rbind(
  cardinalities_footrule_df,
  cardinalities_spearman_df,
  cardinalities_ulam_df
)
cardinalities_df$type <- "cardinalities"
cardinalities_df$message <- "Using exact partition function."

# Then we add importance sampling estimates
load("./data-raw/importance_sampling/estimates.RData")

partition_function_data <- rbind(cardinalities_df, estimates)

# Finally, keep the rows in pfd_old which are not in partition_function_data
partition_function_data <- rbind(
  partition_function_data,
  pfd_old[!interaction(pfd_old$n_items, pfd_old$metric, pfd_old$type) %in%
    interaction(
      partition_function_data$n_items,
      partition_function_data$metric,
      partition_function_data$type
    ), , drop = FALSE]
)

# Finally, save the fit as internal data
usethis::use_data(partition_function_data, internal = TRUE, overwrite = TRUE)
