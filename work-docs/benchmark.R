set.seed(1)
beach_data <- setup_rank_data(
  preferences = beach_preferences
)
system.time({
  model_fit <- compute_mallows(
    data = beach_data,
    compute_options = set_compute_options(save_aug = TRUE))
})
# 0.447 seconds
# 0.33 seconds
