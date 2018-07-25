# Code for generating off-line importance sampling estimates for Spearman and footrule.

# For footrule we have closed form available for n <= 50.
# We therefore compute importance sampling estimates for n > 50.
alpha_vector <- seq(from = 0.01, to = 20, length.out = 50)
n <- 51:56
metric <- "footrule"
nmc <- 1e3

# Compute the estimates in parallel, across n
library(parallel)
cl <- makeCluster(6)

logZ_vector <- parLapply(cl, n,
                         function(n, alpha_vector, metric, nmc) {
                           list(
                             num_items = n,
                             logZ = BayesMallows:::compute_importance_sampling_estimate(alpha_vector, n, metric, nmc)
                             )
                           },
                         alpha_vector = alpha_vector, metric = metric, nmc = nmc)

stopCluster(cl)

# Each element og the list logZ_vector, contains the partition functions for a given n, over alpha_vector
# Perform a polynomial regression to get a fit for each of them
importance_sampling_fit <- lapply(logZ_vector, function(logZ_vector) {
  fit <- lm(logZ ~ poly(alpha, 10),
            data = data.frame(logZ = logZ_vector[['logZ']], alpha = alpha_vector))
  list(
    metric = metric,
    num_items = logZ_vector[['num_items']],
    coefficients = as.numeric(fit$coefficients)
  )
  })


# Finally, save the fit as internal data
devtools::use_data(importance_sampling_fit, internal = TRUE, overwrite = TRUE)
