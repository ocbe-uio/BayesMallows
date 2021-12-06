context("SMC complete rankings: sequence")

#########################
# Generate Dataset
#########################
set.seed(994)

data <- sushi_rankings[1:100,]

# General
n_items <- dim(sushi_rankings)[2]  # Number of items
leap_size <- floor(n_items/5)
metric <- "footrule"

# Generate estimate of Z_n(alpha)
alpha_vector <- seq(from = 0, to = 15, by = 1)
iter <- 1e2
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(method = "importance_sampling",
                                             alpha_vector = alpha_vector,
                                             n_items = n_items, metric = metric,
                                             nmc = iter, degree = degree)


######################################
# BayesMallows Analysis (MCMC)
######################################
nmc <- 20
burnin <- 5
model_fit <- compute_mallows(rankings = data, nmc = nmc, metric = metric, leap_size =leap_size,
                             alpha_prop_sd = 0.15, logz_estimate = logz_estimate)

model_fit$burnin <- burnin

alpha_samples_table <- data.frame(iteration = 1:nmc , value = model_fit$alpha$value)
alpha_samples_table <- alpha_samples_table[(burnin+1):nmc,]

# from observing the plots, this looks like the estimated parameters of the Mallows Model
rho_0 <- c(4,5,2,6,8,3,9,1,7,10)
alpha_0 <- 1.7

# heatplot - there is no burnin!
mcmc_rho_matrix <- matrix(model_fit$rho$value, ncol = n_items, nrow = nmc, byrow=TRUE)

# ###################################################################
# # SMC
# ###################################################################
mcmc_times <- 5
num_new_obs <- 10
Time <- dim(data)[1]/num_new_obs
N <- 100

test <- smc_mallows_new_users_complete(
	R_obs = data, n_items = n_items, metric = metric,
	leap_size = leap_size, N = N, Time = Time,
	logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times,
	alpha_prop_sd = 0.1, lambda = 0.001, alpha_max = 1e6,
	num_new_obs = num_new_obs, verbose = FALSE
)

test_that("Output of smc_mallows_new_users_complete is OK", {
	expect_s3_class(test, "SMCMallows")
	expect_length(test, 2)
	expect_named(test, c("rho_samples", "alpha_samples"))
	expect_equal(dim(test$rho_samples), c(100, 10, 111))
	expect_equal(dim(test$alpha_samples), c(100, 111))
})

# ###############################
# # Analysis
# ###############################

# posterior confidence intervals for rho
rho_temp <- compute_posterior_intervals_rho(
	output = test$rho_samples[,,Time+1], nmc = N, burnin = 0
)

# MAP AND CP consensus ranking estimates
rho_cp <- compute_rho_consensus(
	output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "CP"
)
rho_map <- compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "MAP")

test_that("Output of compute_posterior_intervals_rho is OK", {
	expect_is(rho_temp, "tbl_df")
	expect_length(rho_temp, 7)
	expect_named(
		rho_temp,
		c(
			"item", "parameter", "mean", "median", "conf_level", "hpdi",
			"central_interval"
		)
	)
	expect_equivalent(sapply(rho_temp, length), rep(10, 7))
})

# posterior for alpha
alpha_samples_table <- data.frame(
	iteration = 1:N , value = test$alpha_samples[,Time+1]
)
# posterior confidence intervals
alpha_posterior_intervals <- compute_posterior_intervals_alpha(
	output = test$alpha_samples[,Time+1], nmc = N, burnin = 0
)

test_that("Output of compute_posterior_intervals_alpha is OK", {
	expect_is(alpha_posterior_intervals, "tbl_df")
	expect_length(alpha_posterior_intervals, 6)
	expect_named(
		alpha_posterior_intervals,
		c(
			"parameter", "mean", "median", "conf_level", "hpdi",
			"central_interval"
		)
	)
	expect_equivalent(sapply(alpha_posterior_intervals, length), rep(1, 6))
})

context("SMC complete rankings: breakdown")

test_that("get_mallows_loglik() in smc_mallows_new_users_complete() works", {
	# ======================================================== #
	# Setup                                                    #
	# ======================================================== #

	# Basic elements ----------------------------------------- #
	data <- sushi_rankings[1:100, ]
	n_users <- nrow(data)
	n_items <- ncol(sushi_rankings)
	Time <- nrow(data) / num_new_obs
	num_new_obs <- 10
	N <- 100

	# rho_samples and alpha_samples -------------------------- #
	rho_samples <- array(data=0, dim=c(N, n_items, (n_users + Time + 1)))
	for (ii in seq_len(N)){
		rho_samples[ii, , 1] <- sample(seq_len(n_items), n_items, replace=FALSE)
	}
	alpha_samples <- matrix(nrow=N, ncol=(n_items + Time + 1))
	alpha_samples[, 1] <- rexp(N, rate=1)

	# logz_estimate ------------------------------------------ #
	alpha_vector <- seq(from = 0, to = 15, by = 1)
	iter <- 3e2
	degree <- 10
	logz_estimate <- estimate_partition_function(
		method="importance_sampling", alpha_vector=alpha_vector,
		n_items=n_items, metric=metric, nmc=iter, degree=degree
	)

	num_obs <- 0
	out_loglik <- vector(mode="numeric", length=Time)
	for (tt in seq_len(Time)) {
		num_obs <- num_obs + num_new_obs
		new_observed_rankings <- data[(num_obs - num_new_obs + 1):num_obs, ]
		rho_samples[, , tt + 1] <- rho_samples[, , tt]
		alpha_samples[, tt + 1] <- alpha_samples[, tt]
		alpha_samples_ii <- alpha_samples[ii, tt + 1]
		rho_samples_ii <- rho_samples[ii, , tt + 1]
		for (ii in seq_len(N)) {
			log_z_alpha <- BayesMallows:::get_partition_function(
				n_items, alpha_samples_ii, NULL, logz_estimate, metric
			)
			log_likelihood <- get_mallows_loglik(
				alpha_samples_ii, t(rho_samples_ii), n_items,
				new_observed_rankings, metric
			)
		}
		out_loglik[tt] <- log_likelihood
	}

	# ======================================================== #
	# Test                                                     #
	# ======================================================== #
	tolerance <- 0.1
	expect_gt(max(out_loglik), mean(out_loglik) * (1 + tolerance))
	expect_lt(min(out_loglik), mean(out_loglik) * (1 - tolerance))
})
