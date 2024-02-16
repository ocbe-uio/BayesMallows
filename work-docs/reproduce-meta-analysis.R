library(BayesMallows)
library(RankAggreg)
library(parallel)
data("geneLists")
set.seed(2233)

genes <- sort(unique(as.character(geneLists)))
n_items <- length(genes)
alpha_vector <- seq(from = 0, to = 5, by = .1)


# In JMLR the partition function was estimate for alpha between 0 and 40,
# but since the posterior mean was 0.56, the mode 0.43 and the 95% HPDI was
# (0.04, 1.29), I estimated it instead between 0 and 5.

iterations <- c(1e2, 1e3, 1e4)
estimates <- vector("list", 3)
fits <- vector("list", 3)

cl <- makeCluster(12)

for(i in seq_along(iterations)) {
  fit <- estimate_partition_function(
    method = "importance_sampling",
    alpha_vector = alpha_vector,
    n_items = n_items,
    metric = "footrule",
    n_iterations = iterations[[i]],
    cl = cl
  )
  fits[[i]] <- fit
  estimates[[i]] <-
    vapply(alpha_vector, function(a) sum(a^fit[, 1] * fit[, 2]), numeric(1))
}

stopCluster(cl)

plot(x = alpha_vector,
     y = estimates[[1]],
     type = "l", xlab = expression(alpha), ylab = expression(Z(alpha)))

for(i in seq(from = 2, to = length(iterations))) {
  lines(x = alpha_vector, y = estimates[[i]], col = i)
}

# Create rank data
dat <- matrix(nrow = nrow(geneLists), ncol = n_items)
colnames(dat) <- genes

for(i in seq_len(nrow(dat))) {
  ordering <- match(geneLists[i, ], genes)
  dat[i, ordering] <- seq_along(ordering)
}

cl <- makeCluster(10)
mod <- compute_mallows(
  data = setup_rank_data(rankings = dat),
  model_options = set_model_options(metric = "footrule"),
  compute_options = set_compute_options(
    nmc = 1020000, burnin = 20000, alpha_prop_sd = .5, leap_size = 15,
    alpha_jump = 100, aug_thinning = 100, rho_thinning = 100),
  priors = set_priors(lambda = .1),
  pfun_estimate = fits[[3]],
  cl = cl
)
stopCluster(cl)

assess_convergence(mod)
assess_convergence(mod, parameter = "rho", items = c(34, 4, 50))

plot(mod)
compute_posterior_intervals(mod)

compute_consensus(mod, type = "CP")
heat_plot(mod)
