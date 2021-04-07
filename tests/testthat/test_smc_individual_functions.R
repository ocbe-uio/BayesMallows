context("Testing SMC functions individually")


###############
# test script
###############
# start a new session of R before running this script!

# set.seed(101)
#
# rho = c(1,2,3,4,5,6)
# alpha = 2
# metric = "footrule"
# n_items= 6
#
# get_mallows_loglik(alpha = alpha, rho = rho,  n_items = length(rho), rankings = rho , metric = metric)
# # return 0 because you are comparing the consensus ranking with itself
# # if you change alpha or metric, then the result shall remain as 0
#
# rankings =   sample_mallows(rho0 = rho, alpha0 = alpha, n_samples = 10,
#                             burnin = 1000, thinning = 500)
#
# # depending on your seed, you will get a different collection of rankings in R and C++
# get_mallows_loglik(alpha = alpha, rho = rho,  n_items = n_items, rankings = rankings , metric = metric)
# # -22.6667
test_that("get_mallows_loglik() works as expected", {})


###############
# test script
###############
# # start a new session of R before running this script!
#
# # This functions uses get_mallows_log_lik and leap_and_shift_probs so if the checks match in those worker functions
# # then it is very likely that this function will return the correct outputs.
#
# set.seed(101)
#
# rho = c(1,2,3,4,5,6)
# alpha = 2
# metric = "footrule"
# n_items= 6
#
# rankings =   sample_mallows(rho0 = rho, alpha0 = alpha, n_samples = 10,
#                             burnin = 1000, thinning = 500)
#
# # you can confirm the print statements inside the metropolis_hastings_rho match get_mallows_loglik and leap_and_shift_probs
# test_1 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 1)
# print(test_1)
# # 1 2 3 5 4 6
# # if rho != rho_prime, then it should have a ulam distance of 1
# # if rho == rho_prime, then it should have ulam distance of 0
# BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
# # should be 1
#
#
# test_2 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 2)
# print(test_2)
# # 1 2 3 4 5 6
#
# BayesMallows:::get_rank_distance(rho, test_2, metric= "ulam")
# # should be 0
#
#
# test_3 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rho, metric = metric, rho = rho, leap_size = 3)
# print(test_3)
# # 1 3 4 2 5 6
#
# BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
# # should be 1
#
#
# # we have a ranking data set containing 10 rankings over 6 items
# test_4 = metropolis_hastings_rho(alpha = alpha, n_items = n_items, rankings = rankings, metric = metric, rho = rho, leap_size = 1)
# print(test_4)
# # 1 2 3 5 4 6
#
# BayesMallows:::get_rank_distance(rho, test_4, metric= "ulam")
# # should be 1
#
#
#
test_that("smc_metropolis_hastings_rho() works as expected", {})


########################
# test script
########################
# # start a new R session before running script!
#
# set.seed(101)
#
# rho = c(1,2,3,4,5,6)
# n_items = length(rho)
#
# # leap_size has a possible range, in the BayesMallows papers they suggest leap_size = floor(n_items/5)
# # but the leap_size can be up to n_items/2. Note that leap_size must be integered valued.
# set.seed(101)
#
# # set.seed() will produce different random results in R and C++
#
# # if leap_size = 1, then forwards_prob = backwards_prob
# test_1 = leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 1)
# print(test_1)
# #$rho_prime
# #[1] 1 2 3 4 5 6
#
# #$forwards_prob
# #[1] 0.1666667
#
# #1$backwards_prob
# #[1] 0.1666667
#
# # if rho != rho_prime, then it should have a ulam distance of 1
# # if rho == rho_prime, then it should have ulam distance of 0
# BayesMallows:::get_rank_distance(rho, test_1$rho_prime, metric= "ulam")
# # should be 0
#
#
# test_2 = leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 2)
# print(test_2)
# #$rho_prime
# #[1] 1 2 3 5 6 4
#
# #$forwards_prob
# #[1] 0.08333333
#
# #$backwards_prob
# #[1] 0.04166667
#
# BayesMallows:::get_rank_distance(rho, test_2$rho_prime, metric= "ulam")
# # should be 1
#
#
# test_3 =leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 3)
# print(test_3)
# #$rho_prime
# #[1] 3 1 2 4 5 6
#
# #$forwards_prob
# #[1] 0.05555556
#
# #$backwards_prob
# #[1] 0.03333333
#
# BayesMallows:::get_rank_distance(rho, test_3$rho_prime, metric= "ulam")
# # should be 1
#
test_that("smc_leap_and_shift_probs() works as expected", {})
