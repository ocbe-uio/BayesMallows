context("Testing computation of partition functions")

# Brute force formula
check_log_zn <- function(n, alpha, metric){
  # Generate all permutations
  perm <- gtools::permutations(n, n)

  # Compute the partition function
  if(metric == "footrule") {
    log(sum(exp(- alpha / n * colSums( abs(t(perm ) - 1:n )))))
  } else if(metric == "spearman") {
    log(sum(exp(- alpha / n * colSums( (t(perm ) - 1:n )^2))))
  } else if(metric == "kendall") {
    log(sum(exp(- alpha / n * apply(perm, 1,
                                    get_rank_distance,
                                    r2 = 1:n, metric = "kendall"))))
  } else if(metric == "cayley") {
    log(sum(exp(- alpha / n * apply(perm, 1,
                                    get_rank_distance,
                                    r2 = 1:n, metric = "cayley"))))
  } else {
    stop("Unknown metric.")
  }
}

# Loop over some n and alpha values
test_that("footrule partition function is correct", {
  for(n in c(1, 2, 3, 5)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        get_partition_function(n, alpha, BayesMallows:::footrule_sequence[[n]], "footrule"),
        check_log_zn(n, alpha, "footrule")
      )
    }
  }})

test_that("Spearman partition function is correct", {
  for(n in c(1, 2, 3)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        get_partition_function(n, alpha, BayesMallows:::spearman_sequence[[n]], "spearman"),
        check_log_zn(n, alpha, "spearman")
      )
    }
  }})


test_that("Kendall partition function is correct", {
  for(n in c(1, 2, 3)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        get_partition_function(n, alpha, 0, "kendall"),
        check_log_zn(n, alpha, "kendall")
      )
    }
  }})

test_that("Cayley partition function is correct", {
  for(n in c(1, 2, 3)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        get_partition_function(n, alpha, 0, "cayley"),
        check_log_zn(n, alpha, "cayley")
      )
    }
  }})
