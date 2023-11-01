bmm <- compute_mallows_mixtures(
  n_clusters = c(1, 4),
  data = setup_rank_data(sushi_rankings),
  compute_options = set_compute_options(nmc = 50, include_wcd = FALSE),
  seed = 432
)

expect_equal(
  bmm[[2]]$cluster_probs$value[20:26],
  c(0.109669232480542, 0.364331214168161, 0.152075264807036, 0.39271690112913,
    0.0908766198956731, 0.380642409337839, 0.134505626724183)
)


m <- compute_mallows(
  setup_rank_data(sushi_rankings),
  model = set_model_options(n_clusters = 5),
  seed = 123,
  compute_options = set_compute_options(nmc = 20))

expect_equal(
  m$cluster_probs$value[30:35],
  c(0.747424659240883, 0.100012805539763, 0.0240360776409296, 0.0197784930292096,
    0.0588434949639344, 0.797329128826163)
)
