context("Testing computation of partition functions")

# First we test it for the footrule distance. Remember this is log(Z_n)

# Brute force formula
check_log_zn <- function(n, alpha, metric){
  # Generate all permutations
  perm <- gtools::permutations(n, n)

  # Compute the partition function
  if(metric == "footrule") {
    log(sum(exp(- alpha / n * colSums( abs(t(perm ) - 1:n )))))
  } else if(metric == "spearman") {
    log(sum(exp(- alpha / n * colSums( (t(perm ) - 1:n )^2))))
  }
}

# Loop over some n and alpha values
test_that("footrule partition function is correct", {
  for(n in c(1, 2, 3, 5)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        check_log_zn(n, alpha, "footrule"),
        get_partition_function(n, alpha, BayesMallows:::footrule_sequence[[n]], "footrule")
      )
    }
  }})

test_that("spearman partition function is correct", {
  for(n in c(1, 2, 3)){
    for(alpha in c(0.001, 0.1, 1)){
      expect_equal(
        check_log_zn(n, alpha, "spearman"),
        get_partition_function(n, alpha, BayesMallows:::spearman_sequence[[n]], "spearman")
      )
    }
  }})
