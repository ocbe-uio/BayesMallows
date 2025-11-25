# Asymptotic Approximation of Partition Function

Compute the asymptotic approximation of the logarithm of the partition
function, using the iteration algorithm of Mukherjee (2016) .

## Usage

``` r
asymptotic_partition_function(
  alpha_vector,
  n_items,
  metric,
  K,
  n_iterations = 1000L,
  tol = 1e-09
)
```

## Arguments

- alpha_vector:

  A numeric vector of alpha values.

- n_items:

  Integer specifying the number of items.

- metric:

  One of `"footrule"` and `"spearman"`.

- K:

  Integer.

- n_iterations:

  Integer specifying the number of iterations.

- tol:

  Stopping criterion for algorithm. The previous matrix is subtracted
  from the updated, and if the maximum absolute relative difference is
  below `tol`, the iteration stops.

## Value

A vector, containing the partition function at each value of alpha.

## References

Mukherjee S (2016). “Estimation in exponential families on
permutations.” *The Annals of Statistics*, **44**(2), 853–875.
[doi:10.1214/15-aos1389](https://doi.org/10.1214/15-aos1389) .
