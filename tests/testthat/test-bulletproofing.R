context("SMC bulletproofing")

limit <- sushi_rankings
n_items <- dim(limit)[2]
n_users <- dim(limit)[1]
metric <- "footrule"

leap_size <- floor(n_items / 5)
rho_init <- c(1, 8, 6, 3, 9, 2, 5, 7, 4, 10)
alpha_init <- 5
alpha_prop_sd <- 0.5
lambda <- 0.1
alpha_max <- 20
cluster_rankings <- sushi_rankings[1, ]

# Compute exact logz
logz_list <- prepare_partition_function(metric = metric, n_items = n_items)

test_that("Functions accept rankings and rho as row or column vectors", {
  mhr_r_r <- metropolis_hastings_rho(
    alpha_init, n_items,
    rankings = as.matrix(cluster_rankings),
    rho = rho_init,
    metric,
    leap_size
  )
  mha_r_r <- metropolis_hastings_alpha(
    alpha_init, n_items,
    rankings = as.matrix(cluster_rankings),
    rho = rho_init,
    logz_estimate = NULL,
    cardinalities = logz_list$cardinalities,
    metric,
    alpha_prop_sd,
    alpha_max,
    lambda
  )
  mhr_rt_r <- metropolis_hastings_rho(
    alpha_init, n_items,
    rankings = t(cluster_rankings),
    rho = rho_init,
    metric,
    leap_size
  )
  mha_rt_r <- metropolis_hastings_alpha(
    alpha_init, n_items,
    rankings = t(cluster_rankings),
    rho = rho_init,
    logz_estimate = NULL,
    cardinalities = logz_list$cardinalities,
    metric,
    alpha_prop_sd,
    alpha_max,
    lambda
  )
  mhr_r_rt <- metropolis_hastings_rho(
    alpha_init, n_items,
    rankings = as.matrix(cluster_rankings),
    rho = t(rho_init),
    metric,
    leap_size
  )
  mha_r_rt <- metropolis_hastings_alpha(
    alpha_init, n_items,
    rankings = as.matrix(cluster_rankings),
    rho = t(rho_init),
    logz_estimate = NULL,
    cardinalities = logz_list$cardinalities,
    metric,
    alpha_prop_sd,
    alpha_max,
    lambda
  )
  mhr_rt_rt <- metropolis_hastings_rho(
    alpha_init, n_items,
    rankings = t(cluster_rankings),
    rho = t(rho_init),
    metric,
    leap_size
  )
  mha_rt_rt <- metropolis_hastings_alpha(
    alpha_init, n_items,
    rankings = t(cluster_rankings),
    rho = t(rho_init),
    logz_estimate = NULL,
    cardinalities = logz_list$cardinalities,
    metric,
    alpha_prop_sd,
    alpha_max,
    lambda
  )
  expect_equal(dim(mhr_r_r), c(10, 1))
  expect_length(mha_r_r, 1)
  expect_equal(dim(mhr_rt_r), c(10, 1))
  expect_length(mha_rt_r, 1)
  expect_equal(dim(mhr_r_rt), c(10, 1))
  expect_length(mha_r_rt, 1)
  expect_equal(dim(mhr_rt_rt), c(10, 1))
  expect_length(mha_rt_rt, 1)
})
