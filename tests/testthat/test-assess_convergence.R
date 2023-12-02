test_that("assess_convergence.BayesMallows works for alpha and rho", {
  mod <- compute_mallows(setup_rank_data(potato_visual),
                         compute_options = set_compute_options(nmc = 50))
  p <- assess_convergence(mod)

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$group, "interaction(chain, cluster)")

  p <- assess_convergence(mod, parameter = "rho", items = 1:4)
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")
})

test_that("assess_convergence.BayesMallows works for Rtilde", {
  mod <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE))

  p <- assess_convergence(
    mod, parameter = "Rtilde", items = 1:4, assessors = 1:4)

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")
})

test_that("assess_convergence.BayesMallows works for cluster_probs", {
  mod <- compute_mallows(
    setup_rank_data(rankings = cluster_data),
    compute_options = set_compute_options(nmc = 50),
    model_options = set_model_options(n_clusters = 3)
    )

  p <- assess_convergence(mod, parameter = "cluster_probs")

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "cluster")
})

test_that("assess_convergence.BayesMallows fails properly", {
  mod <- compute_mallows(setup_rank_data(potato_visual),
                         compute_options = set_compute_options(nmc = 3))
  expect_error(
    assess_convergence(mod, parameter = "Rtilde"),
    "Please rerun")

  expect_error(
    assess_convergence(mod, parameter = "alfa"),
    "'arg' should be one of")


})
