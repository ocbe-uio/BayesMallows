context("SMC pseudolikelihood functions")

################################################################################
# test for get_sample_probabilities
################################################################################
set.seed(101)
rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6

item_ordering <- c(3, 6, 4, 5)
partial_ranking <- c(1, 2, NA, NA, NA, NA)
remaining_set <- c(3, 4, 5, 6)
test_1 <- get_sample_probabilities(
  rho_item_rank = rho[3], alpha, remaining_set, metric, n_items
)
test_2 <- get_sample_probabilities(
  rho_item_rank = rho[4], alpha, remaining_set, metric, n_items
)

test_that("get_sample_probabilities outputs as expected", {
  expect_equivalent(test_1, c(0.38493, 0.27581, 0.19763, 0.14161), tol = 1e-4)
  expect_equivalent(test_2, c(0.24318, 0.33938, 0.24318, 0.17424), tol = 1e-4)
})

################################################################################
# test for calculate_forwards_probability and calculate_bacwards_probability
################################################################################
set.seed(104)
rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6

item_ordering <- c(3, 6, 4, 5) # order with which the ranks will be filled
partial_ranking <- c(1, 2, NA, NA, NA, NA)
remaining_set <- c(3, 4, 5, 6)

test_1_forward <- calculate_forward_probability(
  item_ordering = item_ordering, partial_ranking = partial_ranking,
  remaining_set = remaining_set, rho = rho, alpha = alpha,
  n_items = n_items, metric = metric
)
# Tried all combinations of item_ordering. No dice.

current_ranking <- c(1, 2, 6, 5, 4, 3)

test_1_backward_a <- calculate_backward_probability(
  item_ordering = item_ordering, partial_ranking = partial_ranking,
  current_ranking = current_ranking, remaining_set = remaining_set, rho = rho,
  alpha = alpha, n_items = n_items, metric = metric
)

new_current_ranking <- test_1_forward$aug_ranking

test_1_backward_b <- calculate_backward_probability(
  item_ordering = item_ordering, partial_ranking = partial_ranking,
  current_ranking = new_current_ranking, remaining_set = remaining_set, rho = rho,
  alpha = alpha, n_items = n_items, metric = metric
)

test_that("calculations of forward and backward probabilities", {
  expect_equal(test_1_forward$aug_ranking, matrix(c(1, 2, 3, 6, 5, 4)))
  expect_equal(test_1_forward$forward_prob, 0.03699547)
  expect_equal(test_1_backward_a, 0.01360987)
  expect_equal(test_1_backward_b, test_1_forward$forward_prob)
})

############################################
# tests for M-H_aug_ranking_pseudo
###########################################

set.seed(101)
rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6

test_that("M-H aug ranking pseudo works", {
  R_curr <- c(1, 2, 3, 4, 5, 6)
  R_obs <- c(1, 2, 3, 4, 5, 6)
  test_1 <- metropolis_hastings_aug_ranking_pseudo(
    alpha, rho, n_items, R_obs, R_curr, metric
  )
  expect_equal(test_1, matrix(c(1, 2, 3, 4, 5, 6)))
  expect_equal(all(test_1 == R_curr), TRUE)
  R_obs <- c(1, 2, 3, NA, NA, NA)
  set.seed(6220)
  test_2 <- metropolis_hastings_aug_ranking_pseudo(
    alpha, rho, n_items, R_obs, R_curr, metric
  )
  expect_equal(test_2, matrix(c(1, 2, 3, 5, 6, 4)))
  expect_equal(all(test_2 == R_curr), FALSE)
  R_curr <- c(1, 2, 6, 5, 4, 3)
  R_obs <- c(1, 2, NA, NA, NA, NA)
  set.seed(8286)
  test_3 <- metropolis_hastings_aug_ranking_pseudo(
    alpha, rho, n_items, partial_ranking = R_obs, current_ranking = R_curr,
    metric
  )
  expect_equal(test_3, matrix(c(1, 2, 4, 3, 5, 6)))
  expect_equal(all(test_3 == R_curr), FALSE)
})

#########################################################
### tests relating to the correction_kernel function
#########################################################

set.seed(101)
rho <- c(1, 2, 3, 4, 5, 6)
alpha <- 2
metric <- "footrule"
n_items <- 6

test_that("correction_kernel works", {
  # Three missing ------------------------------------------ #
  R_curr <- c(1, 2, 3, 4, 5, 6)
  R_obs <- c(1, 2, 3, NA, NA, NA)
  test_1 <- correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)
  expect_equal(test_1$ranking, as.matrix(c(1, 2, 3, 5, 4, 6)))
  expect_equivalent(test_1$correction_prob, 0.172, tol = 1e-2)
  expect_equal(all(test_1$ranking == R_curr), FALSE)

  # Two missing -------------------------------------------- #
  R_curr <- c(1, 2, 3, 4, 5, 6)
  R_obs <- c(1, 2, 3, 5, NA, NA)
  test_2 <- correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)
  expect_equal(test_2$ranking, as.matrix(c(1, 2, 3, 5, 4, 6)))
  expect_equal(test_2$correction_prob, 0.5)
  expect_equal(all(test_2$ranking == R_curr), FALSE)

  # No missing --------------------------------------------- #
  R_curr <- c(1, 2, 3, 4, 5, 6)
  R_obs <- c(1, 2, 3, 4, 5, 6)
  test_3 <- correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)
  expect_equal(test_3$ranking, as.matrix(c(1, 2, 3, 4, 5, 6)))
  expect_equal(test_3$correction_prob, 1)
  expect_equal(all(test_3$ranking == R_curr), TRUE)
})
