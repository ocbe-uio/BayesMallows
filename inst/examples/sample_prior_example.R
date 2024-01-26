# We can use a collection of particles from the prior distribution as
# initial values for the sequential Monte Carlo algorithm.
# Here we start by drawing 1000 particles from the priors, using default
# parameters.
prior_samples <- sample_prior(1000, ncol(sushi_rankings))
# Next, we provide the prior samples to update_mallws(), together
# with the first five rows of the sushi dataset
model1 <- update_mallows(
  model = prior_samples,
  new_data = setup_rank_data(sushi_rankings[1:5, ]))
plot(model1)

# We keep adding more data
model2 <- update_mallows(
  model = model1,
  new_data = setup_rank_data(sushi_rankings[6:10, ]))
plot(model2)

model3 <- update_mallows(
  model = model2,
  new_data = setup_rank_data(sushi_rankings[11:15, ]))
plot(model3)
