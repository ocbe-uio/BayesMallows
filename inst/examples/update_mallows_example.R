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
burnin(mod_init) <- 2000

# Next, assume we receive four more observations
data_second_batch <- potato_visual[5:8, ]

# We can now update the model using sequential Monte Carlo
mod_second <- update_mallows(
  model = mod_init,
  new_data = setup_rank_data(rankings = data_second_batch),
  smc_options = set_smc_options(resampler = "systematic")
  )

# This model now has a collection of particles approximating the posterior
# distribution after the first and second batch
# We can use all the posterior summary functions as we do for the model
# based on compute_mallows():
plot(mod_second)
plot(mod_second, parameter = "rho", items = 1:4)
compute_posterior_intervals(mod_second)

# Next, assume we receive the third and final batch of data. We can update
# the model again
data_third_batch <- potato_visual[9:12, ]
mod_final <- update_mallows(
  model = mod_second, new_data = setup_rank_data(rankings = data_third_batch))

# We can plot the same things as before
plot(mod_final)
compute_consensus(mod_final)

# UPDATING A MALLOWS MODEL WITH NEW OR UPDATED PARTIAL RANKINGS
# The sequential Monte Carlo algorithm works for data with missing ranks as
# well. This both includes the case where new users arrive with partial ranks,
# and when previously seen users arrive with more complete data than they had
# previously.
# We illustrate for top-k rankings of the first 10 users in potato_visual
potato_top_10 <- ifelse(potato_visual[1:10, ] > 10, NA_real_,
                        potato_visual[1:10, ])
potato_top_12 <- ifelse(potato_visual[1:10, ] > 12, NA_real_,
                        potato_visual[1:10, ])
potato_top_14 <- ifelse(potato_visual[1:10, ] > 14, NA_real_,
                        potato_visual[1:10, ])

# We need the rownames as user IDs
(user_ids <- 1:10)

# First, users provide top-10 rankings
mod_init <- compute_mallows(
  data = setup_rank_data(rankings = potato_top_10, user_ids = user_ids),
  compute_options = set_compute_options(nmc = 10000))

# Convergence seems fine. We set the burnin to 2000.
assess_convergence(mod_init)
burnin(mod_init) <- 2000

# Next assume the users update their rankings, so we have top-12 instead.
mod1 <- update_mallows(
  model = mod_init,
  new_data = setup_rank_data(rankings = potato_top_12, user_ids = user_ids),
  smc_options = set_smc_options(resampler = "stratified")
)

plot(mod1)

# Then, assume we get even more data, this time top-14 rankings:
mod2 <- update_mallows(
  model = mod1,
  new_data = setup_rank_data(rankings = potato_top_14, user_ids = user_ids)
)

plot(mod2)

# Finally, assume a set of new users arrive, who have complete rankings.
potato_new <- potato_visual[11:12, ]
# We need to update the user IDs, to show that these users are different
(user_ids <- 11:12)

mod_final <- update_mallows(
  model = mod2,
  new_data = setup_rank_data(rankings = potato_new, user_ids = user_ids)
)

plot(mod_final)

# We can also update models with pairwise preferences
# We here start by running MCMC on the first 20 assessors of the beach data
# A realistic application should run a larger number of iterations than we
# do in this example.
set.seed(3)
dat <- subset(beach_preferences, assessor <= 20)
mod <- compute_mallows(
  data = setup_rank_data(
    preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 3000, burnin = 1000)
)

# Next we provide assessors 21 to 24 one at a time.
for(i in 21:24){
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(
      preferences = subset(beach_preferences, assessor == i),
      user_ids = i, shuffle_unranked = TRUE),
    smc_options = set_smc_options(latent_sampling_lag = 0)
  )
}

# Compared to running full MCMC, there is a downward bias in the scale
# parameter. This can be alleviated by increasing the number of particles,
# MCMC steps, and the latent sampling lag.
plot(mod)
compute_consensus(mod)
