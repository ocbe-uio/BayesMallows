context("SMC uniform functions")

set.seed(101)
require("BayesMallows")

###############################################
# tests for M-H_aug_ranking function
###############################################

rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6

test_that('MH-aug ranking works', {
	R_curr = c(1,2,3,6,5,4)
	R_obs = c(1,2,3,NA,NA,NA)
	test_1 = metropolis_hastings_aug_ranking(
		current_ranking = R_curr, partial_ranking= R_obs, alpha = alpha,
		rho = rho, n_items = n_items, metric = metric
	)
	expect_equal(test_1, as.matrix(c(1, 2, 3, 5, 6, 4)))
	expect_equal(get_rank_distance(rho, test_1, metric= "ulam"), 1)
	R_curr = rho
	R_obs = c(1,2,3,NA,NA,NA)
	test_2 = metropolis_hastings_aug_ranking(
		current_ranking = R_curr, partial_ranking = R_obs, alpha = alpha,
		rho = rho, n_items = n_items, metric = metric
	)
	expect_equal(test_2, as.matrix(c(1, 2, 3, 4, 5, 6)))
	expect_equal(all(test_2 == rho), TRUE)
	expect_equal(get_rank_distance(rho, test_2, metric= "ulam"), 0)
	R_curr = c(1,2,3,6,5,4)
	R_obs = c(1,2,3,6,5,NA)
	test_3 = metropolis_hastings_aug_ranking(
		current_ranking = R_curr, partial_ranking = R_obs, alpha = alpha,
		rho = rho, n_items = n_items, metric = metric
	)
	expect_equal(test_3, as.matrix(c(1, 2, 3, 6, 5, 4)))
	expect_equal(all(test_3 == R_curr), TRUE)
})

#########################################################
### tests relating to the correction_kernel function
#########################################################

set.seed(101)
n_items= 6

test_that('correction_kernel works', {
	R_curr = c(1,2,3,4,5,6)
	R_obs = c(1,2,3,NA,NA,NA)
	test_4 = correction_kernel(R_curr, R_obs, n_items)
	expect_equal(test_4$ranking, c(1, 2, 3, 4, 5, 6))
	expect_equal(test_4$correction_prob, 1)
	expect_equal(all(test_4$ranking == R_curr), TRUE)
	R_curr = c(1,2,3,4,5,6)
	R_obs = c(1,2,3,5,NA,NA)
	test_5 = correction_kernel(R_curr, R_obs, n_items)
	expect_equal(test_5$ranking, c(1, 2, 3, 5, 4, 6))
	expect_equal(test_5$correction_prob, 0.5)
	expect_equal(all(test_5$ranking == R_curr), FALSE)
	R_curr = c(1,2,3,4,5,6)
	R_obs = c(1,2,3,4,5,6)
	test_6 = correction_kernel(R_curr, R_obs, n_items)
	expect_equal(test_6$ranking, c(1, 2, 3, 4, 5, 6))
	expect_equal(test_6$correction_prob, 1)
	expect_equal(all(test_6$ranking == R_curr), TRUE)
})