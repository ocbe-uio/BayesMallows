set.seed(1)
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 2000, burnin = 500))

heat_plot(model_fit)
heat_plot(model_fit, type = "MAP")
