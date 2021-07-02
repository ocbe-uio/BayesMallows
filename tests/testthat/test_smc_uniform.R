context("SMC uniform functions")
###############
# test script
###############
# # start a new session of R before running this script!
#
# # This functions uses get_mallows_log_lik so if the checks match in the worker function
# # then it is very likely that this function will return the correct outputs.
#

set.seed(101)
require("BayesMallows")

###############################################
# tests for M-H_aug_ranking function
###############################################


rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6


R_curr = c(1,2,3,6,5,4)
R_obs = c(1,2,3,NA,NA,NA)
# test_1 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric) # FIXME: triggering error on get_rank_distance (#98)

# R_curr == rho so we should get rho
R_curr = rho
R_obs = c(1,2,3,NA,NA,NA)
# test_2 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric) # FIXME: triggering error on get_rank_distance (#98)



# R_curr == rho so we should get rho

R_curr = c(1,2,3,6,5,4)
R_obs = c(1,2,3,6,5,NA)
# test_3 = metropolis_hastings_aug_ranking(R_curr = R_curr, R_obs = R_obs, alpha = alpha, rho = rho, n_items = n_items, metric = metric) # FIXME: triggering error on get_rank_distance (#98)

test_that('MH-aug ranking works', {
	# print(test_1)
	# # 1 2 3 4 6 5
	# # if rho != rho_prime, then it should have a ulam distance of 1
	# # if rho == rho_prime, then it should have ulam distance of 0
	# BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
	# # should be 1

	# print(test_2)
	# # 1 2 3 4 5 6

	# all(test_2 == rho)
	# # should be TRUE

	# # if rho != rho_prime, then it should have a ulam distance of 1
	# # if rho == rho_prime, then it should have ulam distance of 0
	# BayesMallows:::get_rank_distance(rho, test_2, metric= "ulam")
	# # should be 0

	# print(test_3)
	# # 1 2 3 6 5 4

	# all(test_3 == R_curr)
	# # should be TRUE
})


#########################################################
### tests relating to the correction_kernel function
#########################################################


set.seed(101)
n_items= 6

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,NA,NA,NA)

test_4 = correction_kernel(R_curr, R_obs, n_items)

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,5,NA,NA)
test_5 = correction_kernel(R_curr, R_obs, n_items)

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,4,5,6)
test_6 = correction_kernel(R_curr, R_obs, n_items)

test_that('correction_kernel works', {
	# print(test_4)
	# #$ranking
	# #[1] 1 2 3 4 5 6
	# #
	# #$correction_prob
	# #[1] 1


	# all(test_4$ranking == R_curr) == TRUE
	# # TRUE

	# print(test_5)
	# #$ranking
	# #[1] 1 2 3 5 4 6
	# #
	# #$correction_prob
	# #[1] 0.5


	# all(test_5$ranking == R_curr) == TRUE
	# #  should be FALSE

	# print(test_6)
	# #$ranking
	# #[1] 1 2 3 4 5 6
	# #
	# #$correction_prob
	# #[1] 1


	# all(test_6$ranking == R_curr) == TRUE
	# #  should be TRUE
})