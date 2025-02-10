get_labs <- function(x) x$labels
if ("get_labs" %in% getNamespaceExports("ggplot2")) {
  get_labs <- ggplot2::get_labs
}

test_that("assess_convergence and plot works for alpha and rho", {
  set.seed(123)
  mod <- compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 50)
  )
  p <- assess_convergence(mod)
  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$group, "interaction(chain, cluster)")
  expect_error(plot(mod), "Please specify the burnin")
  burnin(mod) <- 10
  p <- plot(mod)
  labs <- get_labs(p)
  expect_equal(labs$y, "Posterior density")
  expect_equal(labs$x, expression(alpha))
  expect_error(
    plot(mod, parameter = "alfa"), "'arg' should be one of"
  )
  expect_message(
    p <- plot(mod, parameter = "rho"),
    "Items not provided by user. Picking 5 at random."
  )
  labs <- get_labs(p)
  expect_equal(labs$y, "Posterior probability")
  expect_equal(labs$x, "rank")
  p <- plot(mod, parameter = "rho", items = 1)
  expect_equal(dim(p$data), c(2, 5))
  expect_error(
    plot(mod, parameter = "rho", items = 33),
    "Unknown items."
  )
  expect_error(
    plot(mod, parameter = "rho", items = "A1"),
    "Unknown items."
  )
  p <- plot(mod, parameter = "rho", items = c("P3", "P5"))
  expect_equal(dim(p$data), c(6, 5))

  expect_error(
    p <- assess_convergence(model_fit = mod, parameter = "rho", items = 33:34),
    "numeric items vector must contain indices between 1 and the number of items"
  )
  expect_error(
    p <- assess_convergence(mod, parameter = "rho", items = letters[1:3]),
    "unknown items provided"
  )

  p <- assess_convergence(mod, parameter = "rho", items = 1:4)
  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")

  expect_message(
    p <- assess_convergence(mod, parameter = "rho"),
    "Items not provided by user. Picking 5 at random."
  )
  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")

  mod <- compute_mallows(setup_rank_data(matrix(c(1, 1, 2, 2), ncol = 2)),
    compute_options = set_compute_options(nmc = 5)
  )

  p1 <- get_labs(assess_convergence(mod, parameter = "rho"))
  p2 <- get_labs(assess_convergence(mod, parameter = "rho", items = 1:2))
  p3 <- get_labs(assess_convergence(mod, parameter = "rho", items = 2:1))
  expect_equal(p1, p2)
  expect_equal(p1, p3)
})

test_that("assess_convergence.BayesMallows works for Rtilde", {
  set.seed(123)
  mod <- compute_mallows(
    setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE)
  )

  p <- assess_convergence(
    mod,
    parameter = "Rtilde", items = 1:4, assessors = 1:4
  )

  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")
  expect_error(
    plot(mod, parameter = "Rtilde"),
    "'arg' should be one of"
  )

  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde", items = 1:4),
    "Assessors not provided by user. Picking 5 at random."
  )
  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")

  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde", assessors = 1:4),
    "Items not provided by user. Picking 5 at random."
  )

  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")

  expect_snapshot(p <- assess_convergence(mod, parameter = "Rtilde"))
  labs <- get_labs(p)
  expect_equal(labs$x, "Iteration")
  expect_equal(labs$colour, "item")

  mod <- compute_mallows(
    setup_rank_data(preferences = subset(beach_preferences, assessor <= 3)),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE)
  )

  expect_message(
    p <- assess_convergence(mod, parameter = "Rtilde"),
    "Items not provided by user. Picking 5 at random."
  )

  mod <- compute_mallows(
    setup_rank_data(
      preferences =
        subset(beach_preferences, bottom_item <= 3 & top_item <= 3)
    ),
    compute_options = set_compute_options(nmc = 50, save_aug = TRUE)
  )

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
  set.seed(11)
  mod <- compute_mallows(
    setup_rank_data(rankings = cluster_data),
    compute_options = set_compute_options(nmc = 50, burnin = 10),
    model_options = set_model_options(n_clusters = 3)
  )

  p <- get_labs(assess_convergence(mod, parameter = "rho", items = 1:3))
  expect_equal(p$x, "Iteration")
  expect_equal(p$colour, "item")

  p <- plot(mod, parameter = "cluster_probs")
  expect_equal(dim(p$data), c(120, 4))
  expect_s3_class(p, "ggplot")

  p <- get_labs(assess_convergence(mod, parameter = "cluster_probs"))
  expect_equal(p$x, "Iteration")
  expect_equal(p$colour, "cluster")

  p <- plot(mod, parameter = "cluster_assignment")
  expect_s3_class(p, "ggplot")
  expect_equal(dim(p$data), c(180, 4))
})

test_that("assess_convergence.BayesMallows works for theta", {
  set.seed(123)
  preferences <- data.frame(
    assessor = c(1, 1, 2, 2),
    bottom_item = c(1, 2, 1, 2),
    top_item = c(2, 1, 2, 3)
  )
  mod <- compute_mallows(
    data = setup_rank_data(preferences = preferences),
    model_options = set_model_options(error_model = "bernoulli"),
    compute_options = set_compute_options(nmc = 10, burnin = 2)
  )

  p <- assess_convergence(mod, parameter = "theta")
  expect_equal(get_labs(p)$x, "Iteration")

  p <- plot(mod, parameter = "theta")
  expect_equal(dim(p$data), c(8, 3))
})

test_that("assess_convergence.BayesMallows fails properly", {
  mod <- compute_mallows(setup_rank_data(potato_visual),
    compute_options = set_compute_options(nmc = 3)
  )
  expect_error(
    assess_convergence(mod, parameter = "Rtilde"),
    "Please rerun"
  )

  expect_error(
    assess_convergence(mod, parameter = "alfa"),
    "'arg' should be one of"
  )
})

test_that("assess_convergence.BayesMallowsMixtures works", {
  n_clusters <- seq(from = 1, to = 3)
  models <- compute_mallows_mixtures(
    n_clusters = n_clusters, data = setup_rank_data(cluster_data),
    compute_options = set_compute_options(nmc = 100, include_wcd = TRUE)
  )

  p <- get_labs(assess_convergence(models))
  expect_equal(p$linetype, "Chain")
  expect_equal(p$colour, "Cluster")
  expect_equal(p$x, "Iteration")
  expect_equal(p$group, "interaction(chain, cluster)")

  expect_error(
    assess_convergence(models, parameter = "rho", items = 1:4),
    "'arg' should be one of"
  )

  p <- get_labs(assess_convergence(models, parameter = "cluster_probs"))
  expect_equal(p$x, "Iteration")
  expect_equal(p$colour, "cluster")
})
