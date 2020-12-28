library(dplyr)
library(tidyr)
library(purrr)
context("Testing compute_mallows")

test_that("miscellaneous input validation", {
  expect_error(compute_mallows(nmc = 1000, alpha_prop_sd = 1))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, alpha_jump = 102))
  expect_error(compute_mallows(rankings = potato_visual, lambda = 0))
  expect_error(compute_mallows(rankings = potato_visual, lambda = -10))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, rho_thinning = 200))
  expect_error(compute_mallows(rankings = potato_visual, nmc = 100, aug_thinning = 200))
  expect_error(compute_mallows(rankings = potato_visual, nmc = -100))
})

test_that("rho_init is properly validated",{
  m <- potato_visual
  expect_error(compute_mallows(rankings = m, rho_init = 1:(ncol(m) - 1)))
  expect_error(compute_mallows(rankings = m, rho_init = c(potato_true_ranking[-1], 22)))
  expect_error(compute_mallows(rankings = m, rho_init = c(NA_real_, 2:ncol(m))))
  expect_error(compute_mallows(rankings = m, obs_freq = -1))
  expect_error(compute_mallows(rankings = m, obs_freq = 1))
  expect_error(compute_mallows(rankings = m, obs_freq = 1:11))
}
)

test_that("compute_mallows discovers inconsistent rankings",{
    expect_error(compute_mallows(
      rankings = matrix(c(1, 2, -3,
                          1, 2, 3), nrow = 2, byrow = TRUE)
    ))
  expect_error(compute_mallows(
    rankings = matrix(c(1, 2, 3,
                        1, 2, 2), nrow = 2, byrow = TRUE)
  ))
  }
)


test_that("compute_mallows error model works", {
  preferences <- data.frame(assessor = c(1, 1, 2, 2),
                            bottom_item = c(1, 2, 1, 2),
                            top_item = c(2, 1, 2, 3)
                            )
  expect_error(invisible(capture.output(compute_mallows(preferences = preferences, nmc = 10))))
  expect_s3_class(compute_mallows(preferences = preferences, error_model = "bernoulli", nmc = 10),
                  "BayesMallows")

})

test_that("compute_mallows with missing data works", {
  mat <- potato_visual * ifelse(runif(length(potato_visual)) > 0.8, NA_real_, 1)
  m <- compute_mallows(rankings = mat, nmc = 30)
  expect_gt(sd(m$rho$value), 0)
  expect_gt(sd(m$alpha$value), 0.001)
  expect_s3_class(m, "BayesMallows")

})


test_that("compute_mallows runs with the right distances", {
  for(metric in c("footrule", "spearman", "cayley", "kendall", "ulam", "hamming")){
    expect_s3_class(compute_mallows(potato_visual, metric = metric, nmc = 3), "BayesMallows")
  }

})

test_that("compute_mallows handles integer preferences", {
  m <- beach_preferences %>%
    mutate_all(as.integer) %>%
    compute_mallows(preferences = .)

})

test_that("compute_mallows handles data with lots of missings",{
  R_partial2 <- structure(c(NA, NA, NA, NA, NA, NA, 9, NA, NA, 7, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, 7, 8, 10, NA, NA, 9, NA, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 6, NA, 5, 2,
                            6, 5, 6, 6, 5, 7, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, 3, 4, NA, 3, 3, 3, 3, 4, 5, 3, 3, 3, NA, 3, 3, 4, NA,
                            7, 8, 3, 3, 10, 5, 4, NA, NA, NA, 8, NA, NA, NA, NA, NA, 11,
                            NA, NA, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 4, 2,
                            2, 2, 4, 2, 2, 2, NA, NA, 4, 7, 5, 4, 6, 7, 2, 6, 6, 7, NA, NA,
                            NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9,
                            7, 8, NA, 8, 10, 6, NA, 5, NA, 6, 6, 5, 4, 5, NA, 4, 4, 5, NA,
                            NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, NA, NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                            NA, 9, NA), .Dim = c(12L, 20L))

  m <- compute_mallows(R_partial2)
  expect_s3_class(assess_convergence(m), "gg")

}
          )

test_that("compute_mallows treats obs_freq properly",{
  m1 <- compute_mallows(rankings = potato_visual,
                        obs_freq = rep(1, nrow(potato_visual)), seed = 2233)
  m2 <- compute_mallows(rankings = potato_visual, seed = 2233)
  expect_equal(m1, m2)

  # Test with repeated beach preferences
  obs_freq <- rep(1:4, each = 15)

  # Next, we create a new hypthetical beach_preferences dataframe where each
  # assessor is replicated 1-4 times
  beach_pref_rep <- beach_preferences %>%
    mutate(new_assessor = map(obs_freq[assessor], ~ 1:.x)) %>%
    unnest(cols = new_assessor) %>%
    mutate(assessor = paste(assessor, new_assessor, sep = ",")) %>%
    select(-new_assessor)

  # We generate transitive closure for these preferences
  beach_tc_rep <- generate_transitive_closure(beach_pref_rep)

  # We generate the initial rankings for the repeated and the "unrepeated"
  # data
  set.seed(1223)
  beach_tc <- generate_transitive_closure(beach_preferences)
  beach_rankings <- generate_initial_ranking(beach_tc, n_items = 15)
  beach_rankings_rep <- generate_initial_ranking(beach_tc_rep, n_items = 15)

  model_fit_obs_freq <- compute_mallows(rankings = beach_rankings,
                                       preferences = beach_tc,
                                       obs_freq = obs_freq,
                                       save_aug = TRUE,
                                       nmc = 10, seed = 3344L)

  expect_equal(model_fit_obs_freq$rho$value,
               c(15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15, 3, 12,
                 13, 8, 6, 1, 14, 2, 11, 10, 4, 9, 7, 5, 15, 2, 12, 13, 8, 6,
                 3, 14, 1, 11, 10, 4, 9, 7, 5, 15, 2, 12, 13, 8, 6, 5, 14, 1,
                 11, 10, 3, 9, 7, 4, 15, 2, 12, 13, 8, 6, 5, 14, 1, 11, 10, 3,
                 9, 7, 4, 15, 2, 12, 13, 8, 6, 5, 14, 1, 11, 10, 3, 9, 7, 4, 15,
                 2, 12, 13, 9, 6, 5, 14, 1, 11, 7, 3, 10, 8, 4, 15, 2, 12, 13,
                 8, 6, 5, 14, 1, 11, 7, 3, 10, 9, 4, 15, 2, 11, 13, 8, 6, 5, 14,
                 1, 12, 7, 3, 10, 9, 4, 15, 2, 12, 13, 8, 6, 5, 14, 1, 9, 7, 3,
                 11, 10, 4))

  expect_equal(model_fit_obs_freq$alpha$value,
               c(1, 1, 1, 0.86535647165335, 0.86535647165335, 0.825179803947417,
                 0.751782270910306, 0.751782270910306, 0.751782270910306, 0.751782270910306
               ))

  # Next for the repeated data.
  model_fit_rep <- compute_mallows(rankings = beach_rankings_rep,
                                   preferences = beach_tc_rep,
                                   save_aug = TRUE,
                                   nmc = 10, seed = 3344L)


  expect_equal(model_fit_rep$rho$value,
               c(15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15, 4, 12,
                 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15, 4, 12, 13, 8, 6,
                 2, 14, 3, 11, 10, 1, 9, 7, 5, 14, 4, 12, 13, 8, 6, 2, 15, 3,
                 11, 10, 1, 9, 7, 5, 15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1,
                 9, 7, 5, 15, 4, 12, 13, 8, 6, 2, 14, 3, 11, 10, 1, 9, 7, 5, 15,
                 4, 12, 13, 9, 6, 2, 14, 3, 11, 10, 1, 8, 7, 5, 15, 4, 12, 13,
                 9, 6, 2, 14, 3, 11, 10, 1, 8, 7, 5, 14, 4, 12, 13, 9, 6, 2, 15,
                 3, 11, 10, 1, 8, 7, 5, 14, 4, 12, 13, 9, 6, 2, 15, 3, 11, 10,
                 1, 8, 7, 5))

  expect_equal(model_fit_rep$alpha$value,
               c(1, 1, 0.903043979825357, 0.903043979825357, 0.861700605603424,
                 0.861700605603424, 0.861700605603424, 0.742080897545894, 0.742080897545894,
                 0.742080897545894))
})
