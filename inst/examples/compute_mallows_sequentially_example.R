\dontrun{
# Observe one ranking at each of 12 timepoints
library(ggplot2)
data <- lapply(seq_len(nrow(potato_visual)), function(i) {
  setup_rank_data(potato_visual[i, ], user_ids = i)
})

initial_values <- sample_prior(
  n = 200, n_items = 20,
  priors = set_priors(gamma = 3, lambda = .1))

mod <- compute_mallows_sequentially(
  data = data,
  initial_values = initial_values,
  smc_options = set_smc_options(n_particles = 500, mcmc_steps = 20))

# We can see the acceptance ratio of the move step for each timepoint:
get_acceptance_ratios(mod)

plot_dat <- data.frame(
  n_obs = seq_along(data),
  alpha_mean = apply(mod$alpha_samples, 2, mean),
  alpha_sd = apply(mod$alpha_samples, 2, sd)
)

# Visualize how the dispersion parameter is being learned as more data arrive
ggplot(plot_dat, aes(x = n_obs, y = alpha_mean, ymin = alpha_mean - alpha_sd,
                     ymax = alpha_mean + alpha_sd)) +
  geom_line() +
  geom_ribbon(alpha = .1) +
  ylab(expression(alpha)) +
  xlab("Observations") +
  theme_classic() +
  scale_x_continuous(
    breaks = seq(min(plot_dat$n_obs), max(plot_dat$n_obs), by = 1))

# Visualize the learning of the rank for a given item (item 1 in this example)
plot_dat <- data.frame(
  n_obs = seq_along(data),
  rank_mean = apply(mod$rho_samples[1, , ], 2, mean),
  rank_sd = apply(mod$rho_samples[1, , ], 2, sd)
)

ggplot(plot_dat, aes(x = n_obs, y = rank_mean, ymin = rank_mean - rank_sd,
                     ymax = rank_mean + rank_sd)) +
  geom_line() +
  geom_ribbon(alpha = .1) +
  xlab("Observations") +
  ylab(expression(rho[1])) +
  theme_classic() +
  scale_x_continuous(
    breaks = seq(min(plot_dat$n_obs), max(plot_dat$n_obs), by = 1))
}
