test_that("update_mallows works", {
  set.seed(1)
  data_first_batch <- potato_visual[1:4, ]

  mod_init <- compute_mallows(
    data = setup_rank_data(data_first_batch),
    compute_options = set_compute_options(nmc = 100, burnin = 0)
  )

  data_second_batch <- potato_visual[5:8, ]

  mod_second <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = data_second_batch),
    smc_options = set_smc_options(n_particles = 10)
  )

  expect_error(
    update_mallows(model = mod_init, new_data = data_second_batch),
    "new_data must be an object of class BayesMallowsData"
  )

  pi <- compute_posterior_intervals(mod_second)
  expect_equal(pi$hpdi, "[0.423,0.753]")
  expect_equal(mod_second$alpha$value[[9]], 0.753246865159393)

  data_third_batch <- potato_visual[9:12, ]
  mod_final <- update_mallows(
    model = mod_second, new_data = setup_rank_data(rankings = data_third_batch)
  )

  expect_equal(mod_final$rho$value[169], 18)

  expect_error(
    update_mallows(model = mod_second, new_data = data_third_batch),
    "new_data must be an object of class BayesMallowsData"
  )
})

test_that("update_mallows can start from prior", {
  set.seed(1)
  prior_samples <- sample_prior(100, 20, set_priors(gamma = 2, lambda = .1))
  mod1 <- update_mallows(
    prior_samples,
    new_data = setup_rank_data(potato_visual[1, , drop = FALSE]),
    smc_options = set_smc_options(n_particles = 100)
  )

  expect_error(
    update_mallows(prior_samples, new_data = potato_visual[1, , drop = FALSE]),
    "new_data must be an object of class BayesMallowsData"
  )

  mod2 <- update_mallows(
    mod1,
    new_data = setup_rank_data(potato_visual[2, , drop = FALSE])
  )
  expect_equal(mod2$alpha_samples[[56]], 3.51628380350389)
})

test_that("update_mallows handles estimated partition function", {
  set.seed(199)
  dat <- t(replicate(3, sample(22)))
  fit <- estimate_partition_function(
    method = "asymptotic",
    alpha_vector = seq(from = 0, to = 1, by = .1),
    n_items = 22,
    metric = "spearman",
    n_iterations = 50
  )
  mod <- compute_mallows(
    data = setup_rank_data(dat),
    model_options = set_model_options(metric = "spearman"),
    compute_options = set_compute_options(nmc = 10, burnin = 0),
    pfun_estimate = fit
  )
  expect_equal(mod$pfun_estimate, fit)

  newdat <- t(replicate(3, sample(22)))
  mod1 <- update_mallows(
    model = mod,
    new_data = setup_rank_data(newdat),
    model_options = set_model_options(metric = "spearman"),
    smc_options = set_smc_options(n_particles = 5)
  )
  expect_equal(mod1$pfun_estimate, fit)

  newdat <- t(replicate(3, sample(22)))
  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(newdat)
  )
  expect_equal(mod2$pfun_estimate, fit)
})
