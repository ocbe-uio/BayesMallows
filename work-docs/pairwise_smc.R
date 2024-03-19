library(BayesMallows)
library(tidyverse)
library(patchwork)
#set.seed(1)

dat <- subset(beach_preferences, assessor < 2)

mod_init <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)
# assess_convergence(mod_init)
alpha <- numeric()
rho <- matrix(nrow = 0, ncol = 15)

mod <- mod_init

for(i in 2:60) {
  print(i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(
      preferences = subset(beach_preferences, assessor == i), timepoint = i, user_ids = i),
    smc_options = set_smc_options(
      n_particles = 3000, mcmc_steps = 15, latent_sampling_lag = 2,
      max_topological_sorts = 100)
  )
  alpha <- c(alpha, mean(mod$alpha_samples))
  rho <- rbind(rho, apply(mod$rho_samples, 1, mean))
}

plot(alpha, type = "b")

plot_df <- as.data.frame(rho)
colnames(plot_df) <- paste("Item", 1:15)
plot_df$iteration <- seq_len(nrow(plot_df))

plot_df %>%
  as_tibble() %>%
  pivot_longer(cols = -iteration) %>%
  ggplot(aes(x = iteration, y = value, group = name, color = name)) +
  geom_line()

mod_bmm <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

get_acceptance_ratios(mod_bmm)

plot(mod_bmm) + plot(mod) + plot_layout(ncol = 1)


