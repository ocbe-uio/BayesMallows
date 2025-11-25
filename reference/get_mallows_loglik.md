# Likelihood and log-likelihood evaluation for a Mallows mixture model

Compute either the likelihood or the log-likelihood value of the Mallows
mixture model parameters for a dataset of complete rankings.

## Usage

``` r
get_mallows_loglik(
  rho,
  alpha,
  weights,
  metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
  rankings,
  observation_frequency = NULL,
  log = TRUE
)
```

## Arguments

- rho:

  A matrix of size `n_clusters x n_items` whose rows are permutations of
  the first n_items integers corresponding to the modal rankings of the
  Mallows mixture components.

- alpha:

  A vector of `n_clusters` non-negative scalar specifying the scale
  (precision) parameters of the Mallows mixture components.

- weights:

  A vector of `n_clusters` non-negative scalars specifying the mixture
  weights.

- metric:

  Character string specifying the distance measure to use. Available
  options are `"kendall"`, `"cayley"`, `"hamming"`, `"ulam"`,
  `"footrule"`, and `"spearman"`.

- rankings:

  A matrix with observed rankings in each row.

- observation_frequency:

  A vector of observation frequencies (weights) to apply to each row in
  `rankings`. This can speed up computation if a large number of
  assessors share the same rank pattern. Defaults to `NULL`, which means
  that each row of `rankings` is multiplied by 1. If provided,
  `observation_frequency` must have the same number of elements as there
  are rows in `rankings`, and `rankings` cannot be `NULL`.

- log:

  A logical; if TRUE, the log-likelihood value is returned, otherwise
  its exponential. Default is `TRUE`.

## Value

The likelihood or the log-likelihood value corresponding to one or more
observed complete rankings under the Mallows mixture rank model with
distance specified by the `metric` argument.

## See also

Other rank functions:
[`compute_expected_distance()`](compute_expected_distance.md),
[`compute_observation_frequency()`](compute_observation_frequency.md),
[`compute_rank_distance()`](compute_rank_distance.md),
[`create_ranking()`](create_ranking.md),
[`sample_mallows()`](sample_mallows.md)

## Examples

``` r
# Simulate a sample from a Mallows model with the Kendall distance

n_items <- 5
mydata <- sample_mallows(
  n_samples = 100,
  rho0 = 1:n_items,
  alpha0 = 10,
  metric = "kendall")

# Compute the likelihood and log-likelihood values under the true model...
get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = mydata,
  log = FALSE
  )
#> [1] 1.724993e-68

get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = mydata,
  log = TRUE
  )
#> [1] -156.0306

# or equivalently, by using the frequency distribution
freq_distr <- compute_observation_frequency(mydata)
get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = freq_distr[, 1:n_items],
  observation_frequency = freq_distr[, n_items + 1],
  log = FALSE
  )
#> [1] 1.458226e-47

get_mallows_loglik(
  rho = rbind(1:n_items, 1:n_items),
  alpha = c(10, 10),
  weights = c(0.5, 0.5),
  metric = "kendall",
  rankings = freq_distr[, 1:n_items],
  observation_frequency = freq_distr[, n_items + 1],
  log = TRUE
  )
#> [1] -107.8443
```
