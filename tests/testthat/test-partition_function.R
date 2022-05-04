context("Testing computation of partition functions")

# Brute force formula
check_log_zn <- function(n, alpha, metric) {
  # Generate all permutations
  perm <- permutations(n)


  # Compute the partition function
  if (metric == "footrule") {
    log(sum(exp(-alpha / n * colSums(abs(t(perm) - 1:n)))))
  } else if (metric == "spearman") {
    log(sum(exp(-alpha / n * colSums((t(perm) - 1:n)^2))))
  } else if (metric == "kendall") {
    log(sum(exp(-alpha / n * apply(perm, 1,
      get_rank_distance,
      r2 = 1:n, metric = "kendall"
    ))))
  } else if (metric == "cayley") {
    log(sum(exp(-alpha / n * apply(perm, 1,
      get_rank_distance,
      r2 = 1:n, metric = "cayley"
    ))))
  } else if (metric == "hamming") {
    log(sum(exp(-alpha / n * apply(perm, 1,
      get_rank_distance,
      r2 = 1:n, metric = "hamming"
    ))))
  } else if (metric == "ulam") {
    log(sum(unlist(lapply(seq(0, n - 1, by = 1), function(x) {
      PerMallows::count.perms(perm.length = n, dist.value = x, dist.name = "ulam") * exp(-alpha / n * x)
    }))))
  } else {
    stop("Unknown metric.")
  }
}

# Loop over some n and alpha values
test_that("footrule partition function is correct", {
  footrule_sequence <- subset(
    BayesMallows:::partition_function_data,
    metric == "footrule" & type == "cardinalities"
  )$values

  for (n in c(1, 2, 3, 5)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(
          n = n, alpha = alpha,
          cardinalities = footrule_sequence[[n]], metric = "footrule"
        ),
        check_log_zn(n, alpha, "footrule")
      )
    }
  }
})

test_that("Spearman partition function is correct", {
  spearman_sequence <- subset(
    BayesMallows:::partition_function_data,
    metric == "spearman" & type == "cardinalities"
  )$values

  for (n in c(1, 2, 3)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(
          n = n, alpha = alpha,
          cardinalities = spearman_sequence[[n]], metric = "spearman"
        ),
        check_log_zn(n, alpha, "spearman")
      )
    }
  }
})


test_that("Kendall partition function is correct", {
  for (n in c(1, 2, 3)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(n = n, alpha = alpha, metric = "kendall"),
        check_log_zn(n, alpha, "kendall")
      )
    }
  }
})

test_that("Cayley partition function is correct", {
  for (n in c(1, 2, 3)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(n = n, alpha = alpha, metric = "cayley"),
        check_log_zn(n, alpha, "cayley")
      )
    }
  }
})


test_that("Hamming partition function is correct", {
  for (n in c(1, 2, 3)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(n = n, alpha = alpha, metric = "hamming"),
        check_log_zn(n, alpha, "hamming")
      )
    }
  }
})

test_that("Ulam partition function is correct", {
  ulam_sequence <- subset(
    BayesMallows:::partition_function_data,
    metric == "ulam" & type == "cardinalities"
  )$values
  for (n in c(1, 2, 3)) {
    for (alpha in c(0.001, 0.1, 1)) {
      expect_equal(
        get_partition_function(
          n = n, alpha = alpha,
          cardinalities = ulam_sequence[[n]],
          metric = "ulam"
        ),
        check_log_zn(n, alpha, "ulam")
      )
    }
  }
})


test_that("partition function data is sane", {
  pfd <- aggregate(
    list(n = seq_len(nrow(BayesMallows:::partition_function_data))),
    BayesMallows:::partition_function_data[, c("n_items", "metric", "type")],
    FUN = length)

  expect_equal(
    nrow(subset(pfd, n > 1)),
    0
  )
})

test_that("estimate_partition_function runs in parallel", {
  alpha_vector <- seq(from = 0, to = 10, length.out = 10)
  n_items <- 5
  metric <- "spearman"
  degree <- 2

  fit <- estimate_partition_function(
    method = "importance_sampling",
    alpha_vector = alpha_vector,
    n_items = n_items, metric = metric,
    nmc = 20, degree = degree
  )


  lapply(1:2, function(x) {
    cl <- parallel::makeCluster(x)
    fit <- estimate_partition_function(
      method = "importance_sampling",
      alpha_vector = alpha_vector,
      n_items = n_items, metric = metric,
      nmc = 20, degree = degree, cl = cl
    )
    parallel::stopCluster(cl)
    expect_equal(length(fit), 3)
  })
})
