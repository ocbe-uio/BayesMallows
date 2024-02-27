# Observe one ranking at each of 12 timepoints
data <- lapply(seq_len(nrow(potato_visual))[1:2], function(i) {
  setup_rank_data(potato_visual[i, ])
})

initial_values <- sample_prior(200, 20)
smc_options <- set_smc_options(n_particles = 4)
