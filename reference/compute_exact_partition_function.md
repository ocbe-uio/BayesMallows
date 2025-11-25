# Compute exact partition function

For Cayley, Hamming, and Kendall distances, computationally tractable
functions are available for the exact partition function.

## Usage

``` r
compute_exact_partition_function(
  alpha,
  n_items,
  metric = c("cayley", "hamming", "kendall")
)
```

## Arguments

- alpha:

  Dispersion parameter.

- n_items:

  Number of items.

- metric:

  Distance function, one of "cayley", "hamming", or "kendall".

## Value

The logarithm of the partition function.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## See also

Other partition function:
[`estimate_partition_function()`](estimate_partition_function.md),
[`get_cardinalities()`](get_cardinalities.md)

## Examples

``` r
compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "cayley"
)
#> [1] 85.60544

compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "hamming"
)
#> [1] 85.286

compute_exact_partition_function(
  alpha = 3.4, n_items = 34, metric = "kendall"
)
#> [1] 65.91866
```
