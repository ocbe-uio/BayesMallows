# Set options for Bayesian Mallows model

Specify various model options for the Bayesian Mallows model.

## Usage

``` r
set_model_options(
  metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam"),
  n_clusters = 1,
  error_model = c("none", "bernoulli")
)
```

## Arguments

- metric:

  A character string specifying the distance metric to use in the
  Bayesian Mallows Model. Available options are `"footrule"`,
  `"spearman"`, `"cayley"`, `"hamming"`, `"kendall"`, and `"ulam"`. The
  distance given by `metric` is also used to compute within-cluster
  distances, when `include_wcd = TRUE`.

- n_clusters:

  Integer specifying the number of clusters, i.e., the number of mixture
  components to use. Defaults to `1L`, which means no clustering is
  performed. See
  [`compute_mallows_mixtures()`](compute_mallows_mixtures.md) for a
  convenience function for computing several models with varying numbers
  of mixtures.

- error_model:

  Character string specifying which model to use for inconsistent
  rankings. Defaults to `"none"`, which means that inconsistent rankings
  are not allowed. At the moment, the only available other option is
  `"bernoulli"`, which means that the Bernoulli error model is used. See
  Crispino et al. (2019) for a definition of the Bernoulli model.

## Value

An object of class `"BayesMallowsModelOptions"`, to be provided in the
`model_options` argument to [`compute_mallows()`](compute_mallows.md),
[`compute_mallows_mixtures()`](compute_mallows_mixtures.md), or
[`update_mallows()`](update_mallows.md).

## References

Crispino M, Arjas E, Vitelli V, Barrett N, Frigessi A (2019). “A
Bayesian Mallows approach to nontransitive pair comparison data: How
human are sounds?” *The Annals of Applied Statistics*, **13**(1),
492–519. [doi:10.1214/18-aoas1203](https://doi.org/10.1214/18-aoas1203)
.

## See also

Other preprocessing:
[`get_transitive_closure()`](get_transitive_closure.md),
[`set_compute_options()`](set_compute_options.md),
[`set_initial_values()`](set_initial_values.md),
[`set_priors()`](set_priors.md),
[`set_progress_report()`](set_progress_report.md),
[`set_smc_options()`](set_smc_options.md),
[`setup_rank_data()`](setup_rank_data.md)
