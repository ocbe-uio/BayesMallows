library(BayesMallows)
packageVersion("BayesMallows")
library(tidyverse)
set.seed(1)
leap_size <- 5
example_dataset <- potato_visual
metric <- "footrule"
n_users <- 12
n_items <- 20
test_dataset <- array(0, c(n_users, n_items, (n_items / 2 + 1)))
test_dataset[, , (n_items / 2 + 1)] <- potato_visual
tt <- 0
for (ii in (n_items - 1):(n_items / 2)) {
  tt <- tt + 1

  # set n_users line with one more NA
  example_dataset[example_dataset > ii] <- NA

  # set as new time stamp
  test_dataset[, , ((n_items / 2 + 1) - tt)] <- example_dataset
}

all.equal(as.numeric(test_dataset[,, 11]),
          as.numeric(potato_visual))

alpha_vector <- seq(from = 0, to = 15, by = 0.1)
logz_estimate <- estimate_partition_function(
  method = "importance_sampling",
  alpha_vector = alpha_vector,
  n_items = n_items, metric = "footrule",
  nmc = 1e4, degree = 10
)

Time <- dim(test_dataset)[3]
N <- 5000
aug_method <- "pseudolikelihood"
metric <- "footrule"
smc_test_updated_partial <- smc_mallows_new_item_rank(
  n_items = n_items,
  R_obs = test_dataset,
  metric = metric,
  leap_size = leap_size, N = N,
  Time = Time,
  logz_estimate = logz_estimate,
  mcmc_kernel_app = 20,
  alpha_prop_sd = 0.5,
  lambda = 0.15,
  alpha_max = 1e6,
  aug_method = aug_method
)

bmm_ref <- compute_mallows(rankings = potato_visual, nmc = 2000)
bmm_ref$burnin <- 1000

updated_partial_test_sample_rho <- smc_test_updated_partial$rho_samples[, , Time]
updated_partial_test_sample_alpha <- smc_test_updated_partial$alpha_samples[, Time]
plot_rho_posterior(output = updated_partial_test_sample_rho, nmc = N, burnin = 0, C = 1)
plot_alpha_posterior(output = updated_partial_test_sample_alpha, nmc = N, burnin = 0)

smc_consensus <- compute_rho_consensus(updated_partial_test_sample_rho,
                                       nmc = N, burnin = 0, C = 1, type = "CP")



bmm_consensus <- as.numeric(str_extract(compute_consensus(bmm_ref)$item, "[:digit:]+"))
x <- as.numeric(str_extract(smc_consensus$item, "[:digit:]+"))
y <- bmm_consensus
plot(y ~ x)
abline(lm(y ~ x))


mean(bmm_ref$alpha$value[bmm_ref$alpha$iteration > 1000])
mean(smc_test_updated_partial$alpha_samples[, 11])

median(bmm_ref$alpha$value[bmm_ref$alpha$iteration > 1000])
median(smc_test_updated_partial$alpha_samples[, 11])
