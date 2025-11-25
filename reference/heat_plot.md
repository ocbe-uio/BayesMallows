# Heat plot of posterior probabilities

Generates a heat plot with items in their consensus ordering along the
horizontal axis and ranking along the vertical axis. The color denotes
posterior probability.

## Usage

``` r
heat_plot(model_fit, ...)
```

## Arguments

- model_fit:

  An object of type `BayesMallows`, returned from
  [`compute_mallows()`](compute_mallows.md).

- ...:

  Additional arguments passed on to other methods. In particular,
  `type = "CP"` or `type = "MAP"` can be passed on to
  [`compute_consensus()`](compute_consensus.md) to determine the order
  of items along the horizontal axis.

## Value

A ggplot object.

## Details

In models with a single cluster, the items are sorted along the x-axis
according to the consensus ranking. In models with more than one
cluster, the items are sorted alphabetically.

## See also

Other posterior quantities: [`assign_cluster()`](assign_cluster.md),
[`compute_consensus()`](compute_consensus.md),
[`compute_posterior_intervals()`](compute_posterior_intervals.md),
[`get_acceptance_ratios()`](get_acceptance_ratios.md),
[`plot.BayesMallows()`](plot.BayesMallows.md),
[`plot.SMCMallows()`](plot.SMCMallows.md),
[`plot_elbow()`](plot_elbow.md), [`plot_top_k()`](plot_top_k.md),
[`predict_top_k()`](predict_top_k.md),
[`print.BayesMallows()`](print.BayesMallows.md)

## Examples

``` r
set.seed(1)
model_fit <- compute_mallows(
  setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 2000, burnin = 500))

heat_plot(model_fit)

heat_plot(model_fit, type = "MAP")


## Model with three clusters
mod <- compute_mallows(
  data = setup_rank_data(rankings = cluster_data),
  model_options = set_model_options(n_clusters = 3),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
)

heat_plot(mod)

heat_plot(mod, type = "MAP")
```
