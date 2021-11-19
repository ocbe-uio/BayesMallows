context("Testing computation of distance")

# Brute force formula
check_dist <- function(n, fun){
  # Generate all permutations
  perm <- permutations(n)

  apply(perm, 1, fun, r2 = 1:n)
}

# Loop over some n values
test_that("footrule distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "footrule")
      }),
      check_dist(n, fun = function(r1, r2) sum(abs(r1 - r2)))
      )
  }})

# Loop over some n values
test_that("Spearman distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "spearman")
      }),
      check_dist(n, fun = function(r1, r2) sum((r1 - r2)^2))
    )
  }})

# Loop over some n values
test_that("Kendall distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "kendall")
      }),
      check_dist(n, fun = function(r1, r2)
        PerMallows::distance(r1, r2, "kendall"))
    )
  }})


# Loop over some n values
test_that("Cayley distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "cayley")
      }),
      check_dist(n, fun = function(r1, r2)
        PerMallows::distance(r1, r2, "cayley"))
    )
  }})


# Loop over some n values
test_that("Hamming distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "hamming")
      }),
      check_dist(n, fun = function(r1, r2)
        PerMallows::distance(r1, r2, "hamming"))
    )
  }})


# Loop over some n values
test_that("Ulam distance is correct", {
  for(n in c(2, 3, 5)){
    expect_equal(
      check_dist(n, fun = function(r1, r2) {
        get_rank_distance(r1, r2, "ulam")
      }),
      check_dist(n, fun = function(r1, r2)
        PerMallows::distance(r1, r2, "ulam"))
    )
  }})


test_that("Exported rank_distance is correct", {
  # Distance between two vectors of rankings:
  expect_equal(rank_distance(1:5,5:1, metric = "kendall"), 10)
  expect_equal(
    rank_distance(c(2, 4, 3, 6, 1, 7, 5), c(3, 5, 4, 7, 6, 2, 1), metric = "cayley"),
    6)
  expect_equal(
    rank_distance(c(4, 2, 3, 1), c(3, 4, 1, 2), metric = "hamming"),
    4)
  expect_equal(
    rank_distance(c(1, 3, 5, 7, 9, 8, 6, 4, 2), c(1, 2, 3, 4, 9, 8, 7, 6, 5), "ulam"),
    4)
  expect_equal(
    rank_distance(c(8, 7, 1, 2, 6, 5, 3, 4), c(1, 2, 8, 7, 3, 4, 6, 5), "footrule"),
    32)
  expect_equal(
    rank_distance(c(1, 6, 2, 5, 3, 4), c(4, 3, 5, 2, 6, 1), "spearman"),
    54)

  expect_error(rank_distance(c(1, 6, 2, 5, 3, 4), c(4, 3, 5, 2, 6, 1), "spearman", obs_freq = 1:3))

  expect_equal(
    rank_distance(
      potato_visual, potato_true_ranking, "footrule"
    ),
    c(22, 24, 32, 14, 36, 24, 14, 28, 34, 24, 30, 24)
  )
})
