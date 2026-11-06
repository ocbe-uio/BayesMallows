test_that("Bernoulli error model counts each pairwise comparison once", {
  # Regression test for GitHub issue #446. The shape parameters of the
  # truncated Beta full conditional of theta must equal the prior parameters
  # plus the number of pairwise comparisons which are in disagreement
  # (shape1) and in agreement (shape2) with the current latent rankings.
  set.seed(1)
  dat <- setup_rank_data(
    preferences = bernoulli_data[bernoulli_data$assessor <= 3, ]
  )
  priors <- set_priors()
  nmc <- 20

  fit <- run_mcmc(
    data = dat,
    model_options = set_model_options(error_model = "bernoulli"),
    compute_options = set_compute_options(nmc = nmc, save_aug = TRUE),
    priors = priors,
    initial_values = set_initial_values(),
    pfun_values = extract_pfun_values("footrule", dat$n_items, NULL),
    pfun_estimate = NULL,
    progress_report = set_progress_report()
  )

  prefs <- unique(dat$preferences[, c("assessor", "bottom_item", "top_item")])
  n_comparisons <- nrow(prefs)
  # Number of comparisons contradicting the latent rankings, where the
  # rankings matrix has one column per assessor and one row per item.
  count_mistakes <- function(rankings) {
    sum(rankings[cbind(prefs$bottom_item, prefs$assessor)] <
      rankings[cbind(prefs$top_item, prefs$assessor)])
  }

  # The shape parameters at iteration t are computed from the latent rankings
  # available at the end of iteration t - 1.
  for (t in seq(from = 2, to = nmc)) {
    mistakes <- count_mistakes(fit$augmented_data[, , t - 1])
    expect_equal(fit$shape1[[t]], priors$kappa[[1]] + mistakes)
    expect_equal(fit$shape2[[t]], priors$kappa[[2]] + n_comparisons - mistakes)
  }
  expect_gt(min(fit$shape1[-1]), priors$kappa[[1]])
  expect_lt(max(fit$shape1[-1]), priors$kappa[[1]] + n_comparisons)
})

test_that("Bernoulli error model targets the correct posterior", {
  # Small example for which the posterior distribution can be computed
  # exactly, by enumerating all combinations of the modal ranking rho and the
  # latent rankings of the assessors, integrating out theta analytically and
  # alpha numerically. Regression test for GitHub issue #446.
  preferences <- data.frame(
    assessor = c(1, 1, 1, 2, 2, 2, 3, 3),
    bottom_item = c(2, 3, 1, 2, 3, 3, 1, 3),
    top_item = c(1, 2, 3, 1, 2, 1, 2, 2)
  )
  n_items <- 3
  metric <- "footrule"
  gamma <- 1
  lambda <- 1
  kappa <- c(1, 3)

  perms <- matrix(
    c(1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1),
    ncol = n_items, byrow = TRUE
  )
  n_perms <- nrow(perms)
  dist_mat <- sapply(seq_len(n_perms), function(s) {
    compute_rank_distance(perms, perms[s, ], metric = metric)
  })
  logz <- function(alpha) log(sum(exp(-alpha * dist_mat[, 1] / n_items)))

  prefs <- split(preferences, preferences$assessor)
  n_assessors <- length(prefs)
  n_comparisons <- nrow(preferences)
  # mistakes[r, j] is the number of comparisons made by assessor j which
  # contradict the ranking perms[r, ].
  mistakes <- sapply(prefs, function(p) {
    sapply(seq_len(n_perms), function(r) {
      sum(perms[r, p$bottom_item] < perms[r, p$top_item])
    })
  })

  # All combinations of modal ranking (first column) and latent rankings.
  grid <- as.matrix(expand.grid(rep(list(seq_len(n_perms)), n_assessors + 1)))
  rho_ind <- grid[, 1]
  latent <- grid[, -1, drop = FALSE]
  total_distance <- rowSums(sapply(seq_len(n_assessors), function(j) {
    dist_mat[cbind(latent[, j], rho_ind)]
  }))
  total_mistakes <- rowSums(sapply(seq_len(n_assessors), function(j) {
    mistakes[latent[, j], j]
  }))

  # Integral over alpha of alpha^power times the Mallows likelihood times the
  # gamma prior, for a given total distance d between latent rankings and rho.
  alpha_integral <- function(d, power = 0) {
    integrate(function(alpha) {
      sapply(alpha, function(a) {
        exp(power * log(a) - a * d / n_items - n_assessors * logz(a) +
          dgamma(a, shape = gamma, rate = lambda, log = TRUE))
      })
    }, lower = 0, upper = Inf)$value
  }
  d_values <- sort(unique(total_distance))
  A0 <- sapply(d_values, alpha_integral, power = 0)
  A1 <- sapply(d_values, alpha_integral, power = 1)
  d_index <- match(total_distance, d_values)

  # Integrating out theta, which has a truncated Beta full conditional.
  a <- kappa[[1]] + total_mistakes
  b <- kappa[[2]] + n_comparisons - total_mistakes
  log_weights <- log(A0[d_index]) + lbeta(a, b) +
    pbeta(0.5, a, b, log.p = TRUE)
  weights <- exp(log_weights - max(log_weights))
  weights <- weights / sum(weights)

  exact_rho_probs <- as.numeric(tapply(weights, rho_ind, sum))
  exact_theta_mean <- sum(
    weights * a / (a + b) * pbeta(0.5, a + 1, b) / pbeta(0.5, a, b)
  )
  exact_theta_sq <- sum(
    weights * a * (a + 1) / ((a + b) * (a + b + 1)) *
      pbeta(0.5, a + 2, b) / pbeta(0.5, a, b)
  )
  exact_theta_sd <- sqrt(exact_theta_sq - exact_theta_mean^2)
  exact_alpha_mean <- sum(weights * A1[d_index] / A0[d_index])

  nmc <- 100000
  burnin <- 1000
  set.seed(1)
  mod <- compute_mallows(
    data = setup_rank_data(preferences = preferences),
    model_options = set_model_options(
      error_model = "bernoulli", metric = metric
    ),
    compute_options = set_compute_options(
      nmc = nmc, burnin = burnin, alpha_prop_sd = 0.5
    ),
    priors = set_priors(gamma = gamma, lambda = lambda, kappa = kappa)
  )

  rho_samples <- mod$rho_samples[, 1, -seq_len(burnin)]
  rho_index <- match(
    apply(rho_samples, 2, paste, collapse = ""),
    apply(perms, 1, paste, collapse = "")
  )
  mcmc_rho_probs <- tabulate(rho_index, n_perms) / length(rho_index)
  expect_lt(max(abs(mcmc_rho_probs - exact_rho_probs)), 0.02)

  theta <- mod$theta$value[mod$theta$iteration > burnin]
  expect_equal(mean(theta), exact_theta_mean, tolerance = 0.03)
  expect_equal(sd(theta), exact_theta_sd, tolerance = 0.03)

  alpha <- mod$alpha$value[mod$alpha$iteration > burnin]
  expect_equal(mean(alpha), exact_alpha_mean, tolerance = 0.06)
})
