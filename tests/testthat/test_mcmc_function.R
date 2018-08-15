context("Testing MCMC function on potato data")

# Testing with footrule
model_fit <- compute_mallows(potato_weighing, metric = "footrule", nmc = 10000)
test_that(
  "alpha is in a decent range for footrule",
  expect_true(
    mean(model_fit$alpha) > 10 && mean(model_fit$alpha) < 20
  )
)

test_that(
  "acceptance rate is acceptable for footrule",
  expect_true(
    mean(model_fit$alpha_acceptance) > 0 && mean(model_fit$alpha_acceptance < 1)
  )
)

test_that(
  "acceptance rate is acceptable for rho",
  expect_true(
    mean(model_fit$rho_acceptance) > 0 && mean(model_fit$rho_acceptance) < 1
  )
)

test_that(
  "rho is a rank vector",
  expect_true(
    all(apply(model_fit$rho[, sample.int(n = 5000, size = 100), 1], 2, BayesMallows:::validate_permutation))
  )
)
