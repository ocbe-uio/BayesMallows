test_that("assess_convergence.BayesMallows works for alpha and rho", {
  mod <- compute_mallows(setup_rank_data(potato_visual),
                         compute_options = set_compute_options(nmc = 50))
  p <- assess_convergence(mod)

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$group, "interaction(chain, cluster)")

  expect_error(
    p <- assess_convergence(mod, parameter = "rho", items = 33:34),
    "numeric items vector must contain indices between 1 and the number of items"
  )
  expect_error(
    p <- assess_convergence(mod, parameter = "rho", items = letters[1:3]),
    "unknown items provided"
  )

  p <- assess_convergence(mod, parameter = "rho", items = 1:4)
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  expect_message(
    p <- assess_convergence(mod, parameter = "rho"),
    "Items not provided by user. Picking 5 at random."
  )
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  mod <- compute_mallows(setup_rank_data(matrix(c(1, 1, 2, 2), ncol = 2)),
                         compute_options = set_compute_options(nmc = 5))

  p1 <- assess_convergence(mod, parameter = "rho")
  p2 <- assess_convergence(mod, parameter = "rho", items = 1:2)
  p3 <- assess_convergence(mod, parameter = "rho", items = 2:1)
  expect_equal(p1$labels, p2$labels)
  expect_equal(p1$labels, p3$labels)
})

test_that("assess_convergence.BayesMallows works for Rtilde", {
  mod <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE))

  p <- assess_convergence(
    mod, parameter = "Rtilde", items = 1:4, assessors = 1:4)

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde", items = 1:4),
    "Assessors not provided by user. Picking 5 at random."
  )
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")


  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde", assessors = 1:4),
    "Items not provided by user. Picking 5 at random.")

  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  expect_snapshot(p <- assess_convergence(mod, parameter = "Rtilde"))
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  mod <- compute_mallows(
    setup_rank_data(preferences = subset(beach_preferences, assessor <= 3)),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE))

  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde"),
    "Items not provided by user. Picking 5 at random.")

  mod <- compute_mallows(
    setup_rank_data(
      preferences =
        subset(beach_preferences, bottom_item <= 3 & top_item <= 3)),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE))

  expect_snapshot(p <- assess_convergence(mod, parameter = "Rtilde"))

  expect_error(
    assess_convergence(mod, assessors = 100:103, parameter = "Rtilde"),
    "assessors vector must contain numeric indices between 1 and the number of assessors"
  )

  expect_error(
    assess_convergence(mod, parameter = "theta"),
    "Theta not available. Run compute_mallows with error_model = 'bernoulli'."
  )
})

test_that("assess_convergence.BayesMallows works for cluster_probs", {
  mod <- compute_mallows(
    setup_rank_data(rankings = cluster_data),
    compute_options = set_compute_options(nmc = 50),
    model_options = set_model_options(n_clusters = 3)
    )

  p <- assess_convergence(mod, parameter = "rho", items = 1:3)
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "item")

  p <- assess_convergence(mod, parameter = "cluster_probs")
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "cluster")
})

test_that("assess_convergence.BayesMallows works for theta", {
  preferences <- data.frame(
    assessor = c(1, 1, 2, 2),
    bottom_item = c(1, 2, 1, 2),
    top_item = c(2, 1, 2, 3)
  )
  mod <- compute_mallows(
    data = setup_rank_data(preferences = preferences),
    model_options = set_model_options(error_model = "bernoulli"),
    compute_options = set_compute_options(nmc = 10)
  )

  p <- assess_convergence(mod, parameter = "theta")
  expect_equal(p$labels$x, "Iteration")
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

test_that("assess_convergence.BayesMallowsMixtures works", {
  n_clusters <- seq(from = 1, to = 3)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(cluster_data),
    compute_options = set_compute_options(nmc = 100, include_wcd = TRUE))

  p <- assess_convergence(models)
  expect_equal(p$labels$linetype, "Chain")
  expect_equal(p$labels$colour, "Cluster")
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$group, "interaction(chain, cluster)")

  expect_error(
    assess_convergence(models, parameter = "rho", items = 1:4),
    "'arg' should be one of")

  p <- assess_convergence(models, parameter = "cluster_probs")
  expect_equal(p$labels$x, "Iteration")
  expect_equal(p$labels$colour, "cluster")

})

