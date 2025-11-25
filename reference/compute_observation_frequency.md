# Frequency distribution of the ranking sequences

Construct the frequency distribution of the distinct ranking sequences
from the dataset of the individual rankings. This can be of interest in
itself, but also used to speed up computation by providing the
`observation_frequency` argument to
[`compute_mallows()`](compute_mallows.md).

## Usage

``` r
compute_observation_frequency(rankings)
```

## Arguments

- rankings:

  A matrix with the individual rankings in each row.

## Value

Numeric matrix with the distinct rankings in each row and the
corresponding frequencies indicated in the last `(n_items+1)`-th column.

## See also

Other rank functions:
[`compute_expected_distance()`](compute_expected_distance.md),
[`compute_rank_distance()`](compute_rank_distance.md),
[`create_ranking()`](create_ranking.md),
[`get_mallows_loglik()`](get_mallows_loglik.md),
[`sample_mallows()`](sample_mallows.md)

## Examples

``` r
# Create example data. We set the burn-in and thinning very low
# for the sampling to go fast
data0 <- sample_mallows(rho0 = 1:5, alpha = 10, n_samples = 1000,
                        burnin = 10, thinning = 1)
# Find the frequency distribution
compute_observation_frequency(rankings = data0)
#>       [,1] [,2] [,3] [,4] [,5] [,6]
#>  [1,]    1    2    3    4    5  867
#>  [2,]    1    2    3    5    4   44
#>  [3,]    1    2    4    3    5   27
#>  [4,]    1    3    2    4    5   29
#>  [5,]    1    3    4    2    5    5
#>  [6,]    1    4    2    3    5    7
#>  [7,]    1    4    3    2    5   12
#>  [8,]    2    1    3    4    5    7
#>  [9,]    2    1    3    5    4    2

# The function also works when the data have missing values
rankings <- matrix(c(1, 2, 3, 4,
                     1, 2, 4, NA,
                     1, 2, 4, NA,
                     3, 2, 1, 4,
                     NA, NA, 2, 1,
                     NA, NA, 2, 1,
                     NA, NA, 2, 1,
                     2, NA, 1, NA), ncol = 4, byrow = TRUE)

compute_observation_frequency(rankings)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   NA   NA    2    1    3
#> [2,]    1    2    3    4    1
#> [3,]    1    2    4   NA    2
#> [4,]    2   NA    1   NA    1
#> [5,]    3    2    1    4    1
```
