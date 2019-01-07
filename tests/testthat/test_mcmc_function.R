context("Testing MCMC function on potato data")

library(dplyr)
library(tidyr)

# Testing with footrule
set.seed(200)
model_fit <- compute_mallows(potato_weighing, metric = "footrule", nmc = 1000)
mean_alpha <- pull(summarise(slice(model_fit$alpha, 501:1000), mean(value)))

test_that(
  "alpha is in a decent range for footrule",
  expect_true(mean_alpha > 10 && mean_alpha < 20)
)

test_that(
  "acceptance rate is acceptable for footrule",
  expect_true(
    model_fit$alpha_acceptance > 0 && model_fit$alpha_acceptance < 1
    )
)

test_that(
  "acceptance rate is acceptable for rho",
  expect_true(
    model_fit$rho_acceptance > 0 && model_fit$rho_acceptance < 1
  )
)

test_that(
  "rho is a rank vector",
  expect_true(
    model_fit$rho %>%
      sample_n(1000) %>%
      group_by_all() %>%
      count() %>%
      filter(n > 1) %>%
      nrow() == 0
  )
)
