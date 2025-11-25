# Expected value of metrics under a Mallows rank model

Compute the expectation of several metrics under the Mallows rank model.

## Usage

``` r
compute_expected_distance(
  alpha,
  n_items,
  metric = c("footrule", "spearman", "cayley", "hamming", "kendall", "ulam")
)
```

## Arguments

- alpha:

  Non-negative scalar specifying the scale (precision) parameter in the
  Mallows rank model.

- n_items:

  Integer specifying the number of items.

- metric:

  Character string specifying the distance measure to use. Available
  options are `"kendall"`, `"cayley"`, `"hamming"`, `"ulam"`,
  `"footrule"`, and `"spearman"`.

## Value

A scalar providing the expected value of the `metric` under the Mallows
rank model with distance specified by the `metric` argument.

## See also

Other rank functions:
[`compute_observation_frequency()`](compute_observation_frequency.md),
[`compute_rank_distance()`](compute_rank_distance.md),
[`create_ranking()`](create_ranking.md),
[`get_mallows_loglik()`](get_mallows_loglik.md),
[`sample_mallows()`](sample_mallows.md)

## Examples

``` r
compute_expected_distance(1, 5, metric = "kendall")
#> [1] 4.177277
compute_expected_distance(2, 6, metric = "cayley")
#> [1] 3.212053
compute_expected_distance(1.5, 7, metric = "hamming")
#> [1] 5.761023
compute_expected_distance(5, 30, "ulam")
#> [1] 21.33016
compute_expected_distance(3.5, 45, "footrule")
#> [1] 377.5987
compute_expected_distance(4, 10, "spearman")
#> [1] 7.220669
```
