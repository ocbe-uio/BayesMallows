# Get cardinalities for each distance

The partition function for the Mallows model can be defined in a
computationally efficient manner as \$\$Z\_{n}(\alpha) = \sum\_{d\_{n}
\in \mathcal{D}\_{n}} N\_{m,n} e^{-(\alpha/n) d\_{m}}.\$\$ In this
equation, \\\mathcal{D}\_{n}\\ a set containing all possible distances
at the given number of items, and \\d\_{m}\\ is on element of this set.
Finally, \\N\_{m,n}\\ is the number of possible configurations of the
items that give the particular distance. See Irurozki et al. (2016) ,
Vitelli et al. (2018) , and Crispino et al. (2023) for details.

For footrule distance, the cardinalities come from entry A062869 in the
On-Line Encyclopedia of Integer Sequences (OEIS) (Sloane and Inc. 2020)
. For Spearman distance, they come from entry A175929, and for Ulam
distance from entry A126065.

## Usage

``` r
get_cardinalities(n_items, metric = c("footrule", "spearman", "ulam"))
```

## Arguments

- n_items:

  Number of items.

- metric:

  Distance function, one of "footrule", "spearman", or "ulam".

## Value

A dataframe with two columns, `distance` which contains each distance in
the support set at the current number of items, i.e., \\d\_{m}\\, and
`value` which contains the number of values at this particular
distances, i.e., \\N\_{m,n}\\.

## References

Crispino M, Mollica C, Astuti V, Tardella L (2023). “Efficient and
accurate inference for mixtures of Mallows models with Spearman
distance.” *Statistics and Computing*, **33**(5). ISSN 1573-1375,
[doi:10.1007/s11222-023-10266-8](https://doi.org/10.1007/s11222-023-10266-8)
, <http://dx.doi.org/10.1007/s11222-023-10266-8>.  
  
Irurozki E, Calvo B, Lozano JA (2016). “PerMallows: An R Package for
Mallows and Generalized Mallows Models.” *Journal of Statistical
Software*, **71**(12), 1–30.
[doi:10.18637/jss.v071.i12](https://doi.org/10.18637/jss.v071.i12) .  
  
Sloane NJA, Inc. TOF (2020). “The on-line encyclopedia of integer
sequences.” <https://oeis.org/>.  
  
Vitelli V, Sørensen, Crispino M, Arjas E, Frigessi A (2018).
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research*, **18**(1), 1–49.
<https://jmlr.org/papers/v18/15-481.html>.

## See also

Other partition function:
[`compute_exact_partition_function()`](compute_exact_partition_function.md),
[`estimate_partition_function()`](estimate_partition_function.md)

## Examples

``` r
# Extract the cardinalities for four items with footrule distance
n_items <- 4
dat <- get_cardinalities(n_items)
# Compute the partition function at alpha = 2
alpha <- 2
sum(dat$value * exp(-alpha / n_items * dat$distance))
#> [1] 3.572331
#'
# We can confirm that it is correct by enumerating all possible combinations
all <- expand.grid(1:4, 1:4, 1:4, 1:4)
perms <- all[apply(all, 1, function(x) length(unique(x)) == 4), ]
sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum(abs(x - 1:4)))))
#> [1] 3.572331

# We do the same for the Spearman distance
dat <- get_cardinalities(n_items, metric = "spearman")
sum(dat$value * exp(-alpha / n_items * dat$distance))
#> [1] 2.497585
#'
# We can confirm that it is correct by enumerating all possible combinations
sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum((x - 1:4)^2))))
#> [1] 2.497585
```
