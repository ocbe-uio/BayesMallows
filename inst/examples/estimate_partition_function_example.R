# IMPORTANCE SAMPLING
# Let us estimate logZ(alpha) for 20 items with Spearman distance
# We create a grid of alpha values from 0 to 10
alpha_vector <- seq(from = 0, to = 10, by = 0.5)
n_items <- 20
metric <- "spearman"

# We start with 1e3 Monte Carlo samples
fit1 <- estimate_partition_function(
  method = "importance_sampling", alpha_vector = alpha_vector,
  n_items = n_items, metric = metric, n_iterations = 1e3)
# A matrix containing powers of alpha and regression coefficients is returned
fit1
# The approximated partition function can hence be obtained:
estimate1 <-
  vapply(alpha_vector, function(a) sum(a^fit1[, 1] * fit1[, 2]), numeric(1))

# Now let us recompute with 2e3 Monte Carlo samples
fit2 <- estimate_partition_function(
  method = "importance_sampling", alpha_vector = alpha_vector,
  n_items = n_items, metric = metric, n_iterations = 2e3)
estimate2 <-
  vapply(alpha_vector, function(a) sum(a^fit2[, 1] * fit2[, 2]), numeric(1))

# ASYMPTOTIC APPROXIMATION
# We can also compute an estimate using the asymptotic approximation
fit3 <- estimate_partition_function(
  method = "asymptotic", alpha_vector = alpha_vector,
  n_items = n_items, metric = metric, n_iterations = 50)
estimate3 <-
  vapply(alpha_vector, function(a) sum(a^fit3[, 1] * fit3[, 2]), numeric(1))

# We can now plot the estimates side-by-side
plot(alpha_vector, estimate1, type = "l", xlab = expression(alpha),
     ylab = expression(log(Z(alpha))))
lines(alpha_vector, estimate2, col = 2)
lines(alpha_vector, estimate3, col = 3)
legend(x = 7, y = 40, legend = c("IS,1e3", "IS,2e3", "IPFP"),
       col = 1:3, lty = 1)

# We see that the two importance sampling estimates, which are unbiased,
# overlap. The asymptotic approximation seems a bit off. It can be worthwhile
# to try different values of n_iterations and K.

# When we are happy, we can provide the coefficient vector in the
# pfun_estimate argument to compute_mallows
# Say we choose to use the importance sampling estimate with 1e4 Monte Carlo samples:
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  model_options = set_model_options(metric = "spearman"),
  compute_options = set_compute_options(nmc = 200),
  pfun_estimate = fit2)

