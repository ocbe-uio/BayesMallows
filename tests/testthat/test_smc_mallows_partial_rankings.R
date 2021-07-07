context("SMC new users partial rankings")

#########################
# Generate Dataset
#########################
set.seed(994)

# General
n_items <- dim(sushi_rankings)[2] # Number of items
rho_0 <- seq(from = 1, to = n_items, by = 1) # 'true' consensus ranking
alpha_0 <- 2 # fixed/ 'true' scale parameter
leap_size <- floor(n_items / 5)
metric <- "footrule"

# Generate estimate of Z_n(alpha)
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter <- 1e3
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = metric,
  nmc = iter, degree = degree
)

# Make this information partially observed over time
samples <- sushi_rankings[1:100, ]
samples[samples > 5] <- NA

#######################################
# Bayesmallows MCMC Results
#######################################
nmc <- 2000
bayesmallows_mcmc <- compute_mallows(samples, nmc = nmc, leap_size = leap_size, metric = metric, alpha_prop_sd = 0.15)
bayesmallows_mcmc$burnin <- 1000

# choice items to see in trace plot
items <- sample(1:n_items, 5, replace = F)
items <- sort(items)

test_that('BayesMallows MCMC Results are OK', {
    rho_cp <- compute_consensus(model_fit = bayesmallows_mcmc, type = "CP", burnin = bayesmallows_mcmc$burnin)
    rho_map <- compute_consensus(model_fit = bayesmallows_mcmc, type = "MAP", burnin = bayesmallows_mcmc$burnin)
    post_rho <- compute_posterior_intervals(bayesmallows_mcmc, parameter = "rho")
    post_alpha <- compute_posterior_intervals(bayesmallows_mcmc, parameter = "alpha")
    expect_equal(dim(post_rho), c(10, 7))
    expect_equal(dim(rho_cp), c(10, 3))
    expect_equal(dim(rho_map), c(10, 3))
    expect_equal(dim(post_alpha), c(1, 6))
})

###############################
# SMC Analysis (alpha unknown)
###############################

N <- 100
mcmc_times <- 5
num_new_obs <- 5
Time <- dim(samples)[1] / num_new_obs

alpha_prop_sd <- 0.5
lambda <- 0.15
alpha_max <- 1e6
aug_method <- "random"

test <- smc_mallows_new_users_partial(
  R_obs = samples, n_items = n_items, metric = metric,
  leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times,
  num_new_obs = num_new_obs, alpha_prop_sd = alpha_prop_sd,
  lambda = lambda, alpha_max = alpha_max, aug_method = aug_method
)

test_that('SMC analysis (alpha unknown) results are OK', {
    rho_cp <- compute_rho_consensus(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1, type = "CP")
    rho_map <- compute_rho_consensus(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1, type = "MAP")
    post_rho <- compute_posterior_intervals_rho(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0)
    post_alpha <- compute_posterior_intervals_alpha(output = test$alpha_samples[, Time + 1], nmc = N, burnin = 0)
    expect_equal(dim(rho_cp), c(10, 3))
    expect_equal(dim(rho_map), c(12, 3))
    expect_equal(dim(post_rho), c(10, 7))
    expect_equal(dim(post_alpha), c(1, 6))
})

###############################
# SMC Analysis (alpha fixed)
###############################

alpha_0 <- 2

test <- smc_mallows_new_users_partial_alpha_fixed(
  R_obs = samples, n_items = n_items, metric = metric,
  leap_size = leap_size, N = N, Time = Time,
  logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times,
  num_new_obs = num_new_obs, aug_method = aug_method, alpha = alpha_0
)

test_that('SMC analysis (alpha fixed) results are OK', {
    rho_cp <- compute_rho_consensus(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1, type = "CP")
    rho_map <- compute_rho_consensus(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0, C = 1, type = "MAP")
    post_rho <- compute_posterior_intervals_rho(output = test$rho_samples[, , Time + 1], nmc = N, burnin = 0)
    expect_equal(dim(rho_cp), c(10, 3))
    expect_equal(dim(rho_map), c(12, 3))
    expect_equal(dim(post_rho), c(10, 7))
})