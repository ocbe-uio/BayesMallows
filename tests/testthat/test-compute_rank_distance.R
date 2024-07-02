test_that("compute_rank_distance works", {
  rankings <- 1:5
  rho <- c(1, 3, 2, 5, 4)
  expect_equal(compute_rank_distance(rankings, rho, metric = "footrule"), 4)
  expect_equal(compute_rank_distance(rankings, rho, metric = "spearman"), 4)
  expect_equal(compute_rank_distance(rankings, rho, metric = "cayley"), 2)
  expect_equal(compute_rank_distance(rankings, rho, metric = "hamming"), 4)
  expect_equal(compute_rank_distance(rankings, rho, metric = "kendall"), 2)
  expect_equal(compute_rank_distance(rankings, rho, metric = "ulam"), 2)

  rankings <- matrix(c(1:5, 5:1), byrow = TRUE, ncol = 5)
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "footrule"), c(4, 12)
  )
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "spearman"), c(4, 36)
  )
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "cayley"), c(2, 4)
  )
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "hamming"), c(4, 5)
  )
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "kendall"), c(2, 8)
  )
  expect_equal(
    compute_rank_distance(rankings, rho, metric = "ulam"), c(2, 3)
  )

  expect_error(
    compute_rank_distance(rankings, rho, metric = "miles"),
    "'arg' should be one of"
  )

  observation_frequency <- c(3, 1)
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "footrule",
      observation_frequency = observation_frequency
    ), c(12, 12)
  )
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "spearman",
      observation_frequency = observation_frequency
    ), c(12, 36)
  )
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "cayley",
      observation_frequency = observation_frequency
    ), c(6, 4)
  )
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "hamming",
      observation_frequency = observation_frequency
    ), c(12, 5)
  )
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "kendall",
      observation_frequency = observation_frequency
    ), c(6, 8)
  )
  expect_equal(
    compute_rank_distance(
      rankings, rho,
      metric = "ulam",
      observation_frequency = observation_frequency
    ), c(6, 3)
  )
  expect_equal(
    compute_rank_distance(
      create_ranking(c(5, 1, 4, 3, 2)),
      create_ranking(c(3, 1, 2, 4, 5)),
      "ulam"
    ), 3
  )
})

test_that("distances are right-invariant", {
  n_items <- 10
  rho1 <- sample(n_items)
  rho2 <- sample(n_items)
  eta <- sample(n_items)

  metrics <- c("footrule", "spearman", "kendall", "cayley", "hamming", "ulam")
  names(expectations <- metrics)

  for (metric in metrics) {
    r1 <- compute_rank_distance(rho1, rho2, metric = metric)
    r2 <- compute_rank_distance(rho1[eta], rho2[eta], metric = metric)
    expect_equal(r1, r2)
  }
})
