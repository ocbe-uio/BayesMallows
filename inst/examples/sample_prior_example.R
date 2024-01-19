prior_samples <- sample_prior(10000, ncol(sushi_rankings))
model1 <- update_mallows(
  model = prior_samples,
  new_data = setup_rank_data(sushi_rankings[1:500, ]),
  smc_options = set_smc_options(mcmc_steps = 20))
plot(model1)

model2 <- update_mallows(
  model = model1,
  new_data = setup_rank_data(sushi_rankings[501:1000, ]))
plot(model2)

model3 <- update_mallows(
  model = model2,
  new_data = setup_rank_data(sushi_rankings[1001:1500, ]))
plot(model3)
