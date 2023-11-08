rankings <- matrix(c(1:3, 1:3, 1, 3, 2, 3, 2, 1), ncol = 3, byrow = TRUE)
model_fit <- compute_mallows(
  data = setup_rank_data(rankings = rankings),
  compute_options = set_compute_options(nmc = 3),
  model = set_model_options(n_clusters = 1)
)
