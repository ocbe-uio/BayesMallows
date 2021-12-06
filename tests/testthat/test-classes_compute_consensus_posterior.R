# This script checks the correct usage of the generic functions created during
# the fix of issue 80.
set.seed(6998768)

context("compute_posterior_interval() classes")

# Typical BayesMallows workflow ================================================

fit_bm <- compute_mallows(potato_visual)
fit_bm$burnin <- 1000
fit_bm_post_alpha <- compute_posterior_intervals(fit_bm, parameter = "alpha")
fit_bm_post_rho <- compute_posterior_intervals(fit_bm, parameter = "rho")

# Typical SMC-Mallows workflow =================================================

n_items <- ncol(sushi_rankings)
metric <- "footrule"
alpha_vector <- seq(from = 0, to = 15, by = 0.1)
iter <- 1e3
degree <- 10
logz_estimate <- estimate_partition_function(
  method = "importance_sampling", alpha_vector = alpha_vector,
  n_items = n_items, metric = metric, nmc = iter, degree = degree
)
data <- sushi_rankings[1:100, ]
leap_size <- floor(n_items / 5)
nmc <- N <- 1000
Time <- 20
fit_smc <- smc_mallows_new_users_complete(
  R_obs = data, n_items = n_items, metric = metric, leap_size = leap_size,
  N = N, Time = Time, logz_estimate = logz_estimate, mcmc_kernel_app = 5,
  num_new_obs = 5, alpha_prop_sd = 0.5, lambda = 0.15, alpha_max = 1e6
)
fit_smc_alpha <- fit_smc$alpha_samples[, Time + 1]
fit_smc_post_alpha <- compute_posterior_intervals_alpha(
  output = fit_smc_alpha, nmc = nmc, burnin = 0, verbose = FALSE
)
fit_smc_rho <- fit_smc$rho_samples[, , Time + 1]
fit_smc_post_rho <- compute_posterior_intervals_rho(
  output = fit_smc_rho, nmc = nmc, burnin = 0,
  verbose = FALSE
)

# Emulating the internal workings of compute_posterior_intervals ===============

# BayesMallows ------------------------------------------- #

fit_bm_alpha <- fit_bm$alpha
fit_bm_alpha <- dplyr::group_by(fit_bm_alpha, .data$cluster)
class(fit_bm_alpha) <- c(
  "posterior_BayesMallows", "grouped_df", "tbl_df", "tbl", "data.frame"
)
fit_bm_post_internal_alpha <- .compute_posterior_intervals(
  fit_bm_alpha, "alpha", .95, 3L
)

fit_bm_rho <- fit_bm$rho
fit_bm_rho <- dplyr::group_by(fit_bm_rho, .data$cluster)
class(fit_bm_rho) <- c(
  "posterior_BayesMallows", "grouped_df", "tbl_df", "tbl", "data.frame"
)
fit_bm_post_internal_rho <- .compute_posterior_intervals(
  fit_bm_rho, "rho", .95, 3L
)

# SMC-Mallows -------------------------------------------- #

fit_smc_alpha <- data.frame(iteration = seq_len(nmc), value = fit_smc_alpha)
fit_smc_alpha$n_clusters <- 1
fit_smc_alpha$cluster <- "Cluster 1"
fit_smc_alpha <- dplyr::group_by(fit_smc_alpha, .data$cluster)
class(fit_smc_alpha) <- c(
  "posterior_SMCMallows", "grouped_df", "tbl_df", "tbl", "data.frame"
)
fit_smc_post_internal_alpha <- .compute_posterior_intervals(
  fit_smc_alpha, "alpha", .95, 3L
)

fit_smc_rho <- smc_processing(fit_smc_rho)
fit_smc_rho$n_clusters <- 1
fit_smc_rho$cluster <- "Cluster 1"
fit_smc_rho <- dplyr::group_by(fit_smc_rho, .data$cluster)
class(fit_smc_rho) <- c(
  "posterior_SMCMallows", "grouped_df", "tbl_df", "tbl", "data.frame"
)
fit_smc_post_internal_rho <- .compute_posterior_intervals(
  fit_smc_alpha, "rho", .95, 3L, discrete = TRUE
)

# Testing classes ==============================================================

test_that("Classes are correctly attributed", {
  expect_s3_class(fit_bm, "BayesMallows")
  expect_s3_class(fit_smc, "SMCMallows")
  expect_s3_class(fit_bm_post_alpha, "data.frame")
  expect_s3_class(fit_bm_post_rho, "data.frame")
  expect_s3_class(fit_smc_post_alpha, "data.frame")
  expect_s3_class(fit_smc_post_rho, "data.frame")
  expect_error(.compute_posterior_intervals(fit_bm_post_alpha))
  expect_error(.compute_posterior_intervals(fit_bm_post_rho))
  expect_error(.compute_posterior_intervals(fit_smc_post_alpha))
  expect_error(.compute_posterior_intervals(fit_smc_post_rho))
  expect_s3_class(fit_bm_post_internal_alpha, "data.frame")
  expect_s3_class(fit_bm_post_internal_rho, "data.frame")
  expect_s3_class(fit_smc_post_internal_alpha, "data.frame")
  expect_s3_class(fit_smc_post_internal_rho, "data.frame")
})

context("compute_consensus() classes")

fit_bm_consensus_cp <- compute_consensus(fit_bm, type = "CP")
fit_bm_consensus_map <- compute_consensus(fit_bm, type = "MAP")
fit_smc_rho <- fit_smc$rho_samples[, , Time + 1]
fit_smc_consensus_cp <- compute_rho_consensus(
  output = fit_smc_rho, nmc = nmc, burnin = 0, C = 1, type = "CP"
)
fit_smc_consensus_map <- compute_rho_consensus(
  output = fit_smc_rho, nmc = nmc, burnin = 0, C = 1, type = "MAP"
)

test_that("Classes are correctly attributed", {
  expect_s3_class(fit_bm_consensus_cp, "data.frame")
  expect_s3_class(fit_bm_consensus_map, "data.frame")
  expect_s3_class(fit_smc_consensus_cp, "data.frame")
  expect_s3_class(fit_smc_consensus_map, "data.frame")
})
