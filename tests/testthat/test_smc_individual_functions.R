context("Testing SMC functions individually")

rho <- c(1,2,3,4,5,6)
alpha <- 2
metric <- "footrule"
n_items <- 6

test_that("get_mallows_loglik() works as expected", {
	set.seed(101)
	loglik <- get_mallows_loglik(
		alpha = alpha, rho = rho,  n_items = length(rho), rankings = rho,
		metric = metric
	)
	expect_equal(loglik, 0)

	rankings <- sample_mallows(
		rho0 = rho, alpha0 = alpha, n_samples = 10,
		burnin = 1000, thinning = 500
	)
	loglik <- get_mallows_loglik(
		alpha = alpha, rho = rho,  n_items = n_items, rankings = rankings,
		metric = metric
	)
	expect_equivalent(loglik, -22.6667, tol=1e-4)
})

test_that("smc_metropolis_hastings_rho() works as expected", {
	set.seed(101)
	# This functions uses get_mallows_log_lik and leap_and_shift_probs so if the checks match in those worker functions then it is very likely that this function will return the correct outputs.
	rankings <- sample_mallows(
		rho0 = rho, alpha0 = alpha, n_samples = 10,
		burnin = 1000, thinning = 500
	)

	# you can confirm the print statements inside the metropolis_hastings_rho match get_mallows_loglik and leap_and_shift_probs
	test_1 <- metropolis_hastings_rho(
		alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
		rho = rho, leap_size = 1
	)
	dist_1 <- BayesMallows:::get_rank_distance(rho, test_1, metric= "ulam")
	expect_equal(test_1, c(1, 2, 3, 5, 4, 6))
	# if rho != rho_prime, then it should have a ulam distance of 1
	# if rho == rho_prime, then it should have ulam distance of 0
	expect_equal(dist_1, 1)

	test_2 <- metropolis_hastings_rho(
		alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
		rho = rho, leap_size = 2
	)
	dist_2 <- BayesMallows:::get_rank_distance(rho, test_2, metric = "ulam")
	expect_equal(test_2, c(1, 2, 3, 4, 5, 6))
	expect_equal(dist_2, 0)

	test_3 <- metropolis_hastings_rho(
		alpha = alpha, n_items = n_items, rankings = rho, metric = metric,
		rho = rho, leap_size = 3
	)
	dist_3 <- BayesMallows:::get_rank_distance(rho, test_3, metric = "ulam")
	expect_equal(test_3, c(1, 3, 4, 2, 5, 6)) # FIXME: #74 failing test
	expect_equal(dist_3, 1) # FIXME: #74 failing test

	# we have a ranking data set containing 10 rankings over 6 items
	test_4 <- metropolis_hastings_rho(
		alpha = alpha, n_items = n_items, rankings = rankings, metric = metric,
		rho = rho, leap_size = 1
	)
	dist_4 <- BayesMallows:::get_rank_distance(rho, test_4, metric = "ulam")
	expect_equal(test_4, c(1, 2, 3, 5, 4, 6))
	expect_equal(dist_4, 1)
})

test_that("smc_leap_and_shift_probs() works as expected", {
	# set.seed() will produce different random results in R and C++
	set.seed(101)
	n_items <- length(rho)

	# leap_size has a possible range, the BayesMallows papers suggest leap_size = floor(n_items/5) but the leap_size can be up to n_items/2. Note that leap_size must be integered valued.

	# if leap_size = 1, then forwards_prob = backwards_prob
	test_1 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 1)
	expect_equal(test_1$rho_prime, c(1, 2, 3, 4, 5, 6))
	expect_equivalent(test_1$forwards_prob, 0.1666667, tol=1e-6)
	expect_equivalent(test_1$backwards_prob, 0.1666667, tol=1e-6)

	# if rho != rho_prime, then it should have a ulam distance of 1
	# if rho == rho_prime, then it should have ulam distance of 0
	dist_1 <- BayesMallows:::get_rank_distance(rho, test_1$rho_prime, metric= "ulam")
	expect_equal(dist_1, 0)

	test_2 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 2)
	expect_equal(test_2$rho_prime, c(1, 2, 3, 5, 6, 4))
	expect_equivalent(test_2$forwards_prob, 0.08333333, tol=1e-6)
	expect_equivalent(test_2$backwards_prob, 0.04166667, tol=1e-6)

	dist_2 <- BayesMallows:::get_rank_distance(
		rho, test_2$rho_prime, metric= "ulam"
	)
	expect_equal(dist_2, 1)

	test_3 <- leap_and_shift_probs(rho = rho, n_items = n_items, leap_size = 3)
	expect_equal(test_2$rho_prime, c(3, 1, 2, 4, 5, 6)) # FIXME: #74 fails
	expect_equivalent(test_2$forwards_prob, 0.05555556, tol=1e-6) # FIXME: #74
	expect_equivalent(test_2$backwards_prob, 0.03333333, tol=1e-6) # FIXME: #74

	dist_3 <- BayesMallows:::get_rank_distance(
		rho, test_3$rho_prime, metric= "ulam"
	)
	expect_equal(dist_3, 1)
})
