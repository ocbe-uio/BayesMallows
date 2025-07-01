set.seed(1)
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 2000, burnin = 500))

heat_plot(model_fit)
heat_plot(model_fit, type = "MAP")

## Model with three clusters
mod <- compute_mallows(
  data = setup_rank_data(rankings = cluster_data),
  model_options = set_model_options(n_clusters = 3),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

heat_plot(mod)
heat_plot(mod, type = "MAP")
