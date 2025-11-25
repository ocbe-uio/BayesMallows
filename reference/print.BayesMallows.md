# Print Method for BayesMallows Objects

The default print method for a `BayesMallows` object.

## Usage

``` r
# S3 method for class 'BayesMallows'
print(x, ...)

# S3 method for class 'BayesMallowsMixtures'
print(x, ...)

# S3 method for class 'SMCMallows'
print(x, ...)
```

## Arguments

- x:

  An object of type `BayesMallows`, returned from
  [`compute_mallows()`](compute_mallows.md).

- ...:

  Other arguments passed to `print` (not used).

## See also

Other posterior quantities: [`assign_cluster()`](assign_cluster.md),
[`compute_consensus()`](compute_consensus.md),
[`compute_posterior_intervals()`](compute_posterior_intervals.md),
[`get_acceptance_ratios()`](get_acceptance_ratios.md),
[`heat_plot()`](heat_plot.md),
[`plot.BayesMallows()`](plot.BayesMallows.md),
[`plot.SMCMallows()`](plot.SMCMallows.md),
[`plot_elbow()`](plot_elbow.md), [`plot_top_k()`](plot_top_k.md),
[`predict_top_k()`](predict_top_k.md)
