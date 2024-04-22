devtools::load_all()
library(ggplot2)

dat <- lapply(1:150, function(i) {
  setup_rank_data(sushi_rankings[i, ], user_ids = i)
})

initial_values <- sample_prior(
  n = 3000, n_items = 10,
  priors = set_priors(gamma = 3, lambda = .1)
)

mod <- compute_mallows_sequentially(
  data = dat,
  initial_values = initial_values,
  smc_options = set_smc_options(n_particles = 500, resampling_threshold = 250)
)

plot_dat <- data.frame(
  n_obs = seq_along(dat),
  alpha_mean = apply(mod$alpha_samples, 2, mean),
  alpha_sd = apply(mod$alpha_samples, 2, sd)
)

ggplot(plot_dat, aes(x = n_obs, y = alpha_mean, ymin = alpha_mean - alpha_sd,
                     ymax = alpha_mean + alpha_sd)) +
  geom_line() +
  geom_ribbon(alpha = .1) +
  ylab(expression(alpha)) +
  xlab("Observations") +
  theme_classic()

plot_dat <- data.frame(
  n_obs = seq_along(dat),
  rank_mean = apply(mod$rho_samples[3, , ], 2, mean),
  rank_sd = apply(mod$rho_samples[3, , ], 2, sd)
)

ggplot(plot_dat, aes(x = n_obs, y = rank_mean, ymin = rank_mean - rank_sd,
                     ymax = rank_mean + rank_sd)) +
  geom_line() +
  geom_ribbon(alpha = .1) +
  xlab("Observations") +
  ylab(expression(rho[3])) +
  theme_classic()
