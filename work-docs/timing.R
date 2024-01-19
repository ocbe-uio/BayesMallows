library(BayesMallows)
set.seed(1)
# UPDATING A MALLOWS MODEL WITH NEW COMPLETE RANKINGS
# Assume we first only observe the first four rankings in the potato_visual
# dataset
data_first_batch <- potato_visual[1:4, ]

# We start by fitting a model using Metropolis-Hastings
mod_init <- compute_mallows(
  data = setup_rank_data(data_first_batch),
  compute_options = set_compute_options(nmc = 10000))

# Convergence seems good after no more than 2000 iterations
assess_convergence(mod_init)
mod_init$burnin <- 2000

# Next, assume we receive four more observations
data_second_batch <- potato_visual[5:8, ]

# We can now update the model using sequential Monte Carlo
microbenchmark::microbenchmark({
  mod_second <- update_mallows(
    model = mod_init, new_data = setup_rank_data(rankings = data_second_batch))

})
# 14.5 ms
# 9.6 ms in parallel
