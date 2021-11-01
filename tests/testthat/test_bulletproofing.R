context("SMC bulletproofing")

limit <- sushi_rankings
n_items <- dim(limit)[2]
n_users <- dim(limit)[1]
metric <- "footrule"

# Generate estimate of Z_n(alpha)
# Estimate the logarithm of the partition function of the Mallows rank model.
# We create a grid of alpha values from 0 to 10
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e1
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(
  "importance_sampling", alpha_vector, n_items, metric, iter, degree
)
leap_size <- floor(n_items / 5)
rho_init <- c(1, 8, 6, 3, 9, 2, 5, 7, 4, 10)
alpha_init <- 5
alpha_prop_sd <- 0.5
lambda <- 0.1
alpha_max <- 20
cluster_rankings <- sushi_rankings[1, ]

test_that("Functions accept rankings and rho as row or column vectors", {
  mhr1 <- metropolis_hastings_rho(
    alpha_init, n_items,
    rankings = t(as.matrix(cluster_rankings)),
    metric,
    rho = rho_init,
    leap_size
  )
  mha1 <- metropolis_hastings_alpha(
    alpha_init, n_items,
    rankings = t(as.matrix(cluster_rankings)),
    metric,
    rho = rho_init,
    logz_estimate,
    alpha_prop_sd,
    lambda,
    alpha_max
  )
  expect_equal(dim(mhr1), c(10, 1))
  expect_length(mha1, 1)
})
