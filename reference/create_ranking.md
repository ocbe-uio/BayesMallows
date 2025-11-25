# Convert between ranking and ordering.

`create_ranking` takes a vector or matrix of ordered items `orderings`
and returns a corresponding vector or matrix of ranked items.
`create_ordering` takes a vector or matrix of rankings `rankings` and
returns a corresponding vector or matrix of ordered items.

## Usage

``` r
create_ranking(orderings)

create_ordering(rankings)
```

## Arguments

- orderings:

  A vector or matrix of ordered items. If a matrix, it should be of size
  N times n, where N is the number of samples and n is the number of
  items.

- rankings:

  A vector or matrix of ranked items. If a matrix, it should be N times
  n, where N is the number of samples and n is the number of items.

## Value

A vector or matrix of rankings. Missing orderings coded as `NA` are
propagated into corresponding missing ranks and vice versa.

## See also

Other rank functions:
[`compute_expected_distance()`](compute_expected_distance.md),
[`compute_observation_frequency()`](compute_observation_frequency.md),
[`compute_rank_distance()`](compute_rank_distance.md),
[`get_mallows_loglik()`](get_mallows_loglik.md),
[`sample_mallows()`](sample_mallows.md)

## Examples

``` r
# A vector of ordered items.
orderings <- c(5, 1, 2, 4, 3)
# Get ranks
rankings <- create_ranking(orderings)
# rankings is c(2, 3, 5, 4, 1)
# Finally we convert it backed to an ordering.
orderings_2 <- create_ordering(rankings)
# Confirm that we get back what we had
all.equal(orderings, orderings_2)
#> [1] TRUE

# Next, we have a matrix with N = 19 samples
# and n = 4 items
set.seed(21)
N <- 10
n <- 4
orderings <- t(replicate(N, sample.int(n)))
# Convert the ordering to ranking
rankings <- create_ranking(orderings)
# Now we try to convert it back to an ordering.
orderings_2 <- create_ordering(rankings)
# Confirm that we get back what we had
all.equal(orderings, orderings_2)
#> [1] TRUE
```
