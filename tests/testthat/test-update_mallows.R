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

  pi <- compute_posterior_intervals(mod_second)
  expect_equal(pi$hpdi, "[0.331,0.771]")
  expect_equal(mod_second$alpha$value[[9]], 0.528698473629568)

  data_third_batch <- potato_visual[9:12, ]
  mod_final <- update_mallows(
    model = mod_second, new_data = setup_rank_data(rankings = data_third_batch)
  )

  expect_equal(mod_final$rho$value[169], 17)

  skip_on_os("mac")
  potato_top_10 <- ifelse(potato_visual[1:10, ] > 10, NA_real_,
    potato_visual[1:10, ]
  )
  potato_top_12 <- ifelse(potato_visual[1:10, ] > 12, NA_real_,
    potato_visual[1:10, ]
  )
  potato_top_14 <- ifelse(potato_visual[1:10, ] > 14, NA_real_,
    potato_visual[1:10, ]
  )

  user_ids <- rownames(potato_visual[1:10, ])

  mod_init <- compute_mallows(
    data = setup_rank_data(rankings = potato_top_10, user_ids = user_ids),
    compute_options = set_compute_options(nmc = 100, burnin = 5)
  )

  mod1 <- update_mallows(
    model = mod_init,
    new_data = setup_rank_data(rankings = potato_top_12, user_ids = user_ids),
    smc_options = set_smc_options(n_particles = 20)
  )

  expect_equal(mod1$alpha$value[[13]], 0.544567207683911)

  mod2 <- update_mallows(
    model = mod1,
    new_data = setup_rank_data(rankings = potato_top_14, user_ids = user_ids)
  )
  expect_equal(mod2$rho$value[[300]], 7)

  potato_new <- potato_visual[11:12, ]
  user_ids <- rownames(potato_new)

  mod_final <- update_mallows(
    model = mod2,
    new_data = setup_rank_data(rankings = potato_new, user_ids = user_ids)
  )
  expect_equal(mod_final$rho$value[[300]], 4)
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
