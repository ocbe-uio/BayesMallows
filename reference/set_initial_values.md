# Set initial values of scale parameter and modal ranking

Set initial values used by the Metropolis-Hastings algorithm.

## Usage

``` r
set_initial_values(rho_init = NULL, alpha_init = 1)
```

## Arguments

- rho_init:

  Numeric vector specifying the initial value of the latent consensus
  ranking \\\rho\\. Defaults to NULL, which means that the initial value
  is set randomly. If `rho_init` is provided when `n_clusters > 1`, each
  mixture component \\\rho\_{c}\\ gets the same initial value.

- alpha_init:

  Numeric value specifying the initial value of the scale parameter
  \\\alpha\\. Defaults to `1`. When `n_clusters > 1`, each mixture
  component \\\alpha\_{c}\\ gets the same initial value. When chains are
  run in parallel, by providing an argument `cl = cl`, then `alpha_init`
  can be a vector of of length `length(cl)`, each element of which
  becomes an initial value for the given chain.

## Value

An object of class `"BayesMallowsInitialValues"`, to be provided to the
`initial_values` argument of [`compute_mallows()`](compute_mallows.md)
or [`compute_mallows_mixtures()`](compute_mallows_mixtures.md).

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_compute_options()`](set_compute_options.md),
[`set_model_options()`](set_model_options.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)
