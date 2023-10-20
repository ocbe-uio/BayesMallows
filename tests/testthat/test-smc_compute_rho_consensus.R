context("Further tests on compute_rho_consensus")
set.seed(1234)

####################################
# generate example for testing
####################################
n_items <- dim(sushi_rankings)[2]
metric <- "footrule"
alpha_vector <- seq(from = 0, to = 15, by = 0.1)
iter <- 1e2
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)

data <- sushi_rankings[1:100, ]
leap_size <- floor(n_items / 5)
N <- 100
Time <- 20

model_fit <- smc_mallows_new_users(
  R_obs = data, "complete", n_items = n_items,
  metric = metric, leap_size = leap_size,
  N = N, Time = Time,
  mcmc_kernel_app = 5,
  num_new_obs = 5,
  alpha_prop_sd = 0.5,
  lambda = 0.15,
  alpha_max = 1e6
)

##################
# Run tests
##################


test1 <- compute_consensus(model_fit, type = "CP")
test2 <- compute_consensus(model_fit, type = "MAP")

test_that("Output of compute_rho_consensus (CP) is OK", {
  expect_is(test1, "data.frame")
  expect_length(test1, 3)
  expect_named(test1, c("ranking", "item", "cumprob"))
  expect_equal(dim(test1), c(n_items, 3))
  expect_equal(sum(table(test1$ranking)), n_items)
  expect_equal(sum(table(test1$item)), n_items)
})

test_that("Output of compute_rho_consensus (MAP) is OK", {
  expect_is(test2, "data.frame")
  expect_length(test2, 3)
  expect_named(test2, c("probability", "item", "map_ranking"))
  expect_equal(dim(test2), c(n_items, 3))
  expect_equal(sum(table(test2$map_ranking)), n_items)
  expect_equal(sum(table(test2$item)), n_items)
})
