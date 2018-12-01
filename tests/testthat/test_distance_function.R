context("Testing computation of distance")

# Brute force formula
check_dist <- function(n, fun){
  # Generate all permutations
  perm <- gtools::permutations(n, n)

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
