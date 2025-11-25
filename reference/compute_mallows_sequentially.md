# Estimate the Bayesian Mallows Model Sequentially

Compute the posterior distributions of the parameters of the Bayesian
Mallows model using sequential Monte Carlo. This is based on the
algorithms developed in Stein (2023) . This function differs from
[`update_mallows()`](update_mallows.md) in that it takes all the data at
once, and uses SMC to fit the model step-by-step. Used in this way, SMC
is an alternative to Metropolis-Hastings, which may work better in some
settings. In addition, it allows visualization of the learning process.

## Usage

``` r
compute_mallows_sequentially(
  data,
  initial_values,
  model_options = set_model_options(),
  smc_options = set_smc_options(),
  compute_options = set_compute_options(),
  priors = set_priors(),
  pfun_estimate = NULL
)
```

## Arguments

- data:

  A list of objects of class "BayesMallowsData" returned from
  [`setup_rank_data()`](setup_rank_data.md). Each list element is
  interpreted as the data belonging to a given timepoint.

- initial_values:

  An object of class "BayesMallowsPriorSamples" returned from
  [`sample_prior()`](sample_prior.md).

- model_options:

  An object of class "BayesMallowsModelOptions" returned from
  [`set_model_options()`](set_model_options.md).

- smc_options:

  An object of class "SMCOptions" returned from
  [`set_smc_options()`](set_smc_options.md).

- compute_options:

  An object of class "BayesMallowsComputeOptions" returned from
  [`set_compute_options()`](set_compute_options.md).

- priors:

  An object of class "BayesMallowsPriors" returned from
  [`set_priors()`](set_priors.md).

- pfun_estimate:

  Object returned from
  [`estimate_partition_function()`](estimate_partition_function.md).
  Defaults to `NULL`, and will only be used for footrule, Spearman, or
  Ulam distances when the cardinalities are not available, cf.
  [`get_cardinalities()`](get_cardinalities.md).

## Value

An object of class BayesMallowsSequential.

## Details

This function is very new, and plotting functions and other tools for
visualizing the posterior distribution do not yet work. See the examples
for some workarounds.

## References

Stein A (2023). *Sequential Inference with the Mallows Model*. Ph.D.
thesis, Lancaster University.

## See also

Other modeling: [`burnin()`](burnin.md), `burnin<-()`,
[`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md),
[`sample_prior()`](sample_prior.md),
[`update_mallows()`](update_mallows.md)

## Examples

``` r
if (FALSE) { # \dontrun{
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
} # }
```
