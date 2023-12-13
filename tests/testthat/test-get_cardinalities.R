test_that("get_cardinalities works", {
  n_items <- 5
  s <- seq(from = 1, to = n_items, by = 1)
  all <- expand.grid(s, s, s, s, s)
  perms <- all[apply(all, 1, function(x) {length(unique(x)) == n_items}),]

  alpha <- 2
  dat <- get_cardinalities(n_items)

  expect_equal(
    sum(dat$value * exp(-alpha / n_items * dat$distance)),
    sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum(abs(x - s)))))
  )
  expect_equal(
    sum(dat$value * exp(-alpha / n_items * dat$distance)),
    sum(exp(-alpha / n_items * compute_rank_distance(as.matrix(perms), s)))
  )

  dat <- get_cardinalities(n_items, metric = "spearman")
  expect_equal(
    sum(dat$value * exp(-alpha / n_items * dat$distance)),
    sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum((x - s)^2))))
  )
  expect_equal(
    sum(dat$value * exp(-alpha / n_items * dat$distance)),
    sum(exp(-alpha / n_items *
              compute_rank_distance(as.matrix(perms), s, metric = "spearman")))
  )

  dat <- get_cardinalities(n_items, metric = "ulam")
  expect_equal(
    sum(dat$value * exp(-alpha / n_items * dat$distance)),
    sum(exp(-alpha / n_items *
              compute_rank_distance(as.matrix(perms), s, metric = "ulam")))
  )
})
