`BayesMallows`: An R Package for Probabilistic Preference Learning with
the Mallows Rank Model
================
Ã˜ystein SÃ¸rensen
2018-08-13

The `BayesMallows` package implements methods for Bayesian preference
learning with the Mallows rank model, as originally described in Vitelli
et al. (2018), and further developed in Asfaw et al. (2016) and Crispino
et al. (2018). This vignette describes the usage of the package,
starting from the complete data cases, through top-\(k\) rankings,
pairwise comparisons, and finally clustering. We refer to the above
mentioned papers, as well as the review Liu et al. (2018) for a thorough
description of the methods. The necessary methods for data
preprocessing, tuning of algorithms, and assessment of the posterior
distributions will be described along the way.

[The package’s GitHub
repository](https://github.uio.no/oyss/BayesMallows) describes how to
install it. We start this vignette by loading it.

``` r
library(BayesMallows)
```

# Completely Ranked Data

## Potato Data

`BayesMallows` comes with example data described in Liu et al. (2018). A
total of 12 assessors were asked to rank 20 potatoes based on their
weight. In the first round, the assessors were only allowed to study the
potatoes visually, while in the second round, the assessors were also
allowed to hold the potatoes in their hands in order to compare them.
The data sets are named `potato_visual` and `potato_weighing`,
respectively. The true ordering of the potatoes’ weights are stored in
the vector `potato_true_ranking`.

The `potato_visual` dataset is shown below. The column names P1, …, P20
represent potatoes, and the row names A1, …, A12 represent assessors.
The `potato_weighing` dataset has a similar
structure.

|     | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 | P10 | P11 | P12 | P13 | P14 | P15 | P16 | P17 | P18 | P19 | P20 |
| --- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| A1  | 10 | 18 | 19 | 15 |  6 | 16 |  4 | 20 |  3 |   5 |  12 |   1 |   2 |   9 |  17 |   8 |   7 |  14 |  13 |  11 |
| A2  | 10 | 18 | 19 | 17 | 11 | 15 |  6 | 20 |  4 |   3 |  13 |   1 |   2 |   7 |  16 |   8 |   5 |  12 |   9 |  14 |
| A3  | 12 | 15 | 18 | 16 | 13 | 11 |  7 | 20 |  6 |   3 |   8 |   2 |   1 |   4 |  19 |   5 |   9 |  14 |  10 |  17 |
| A4  |  9 | 17 | 19 | 16 | 10 | 15 |  5 | 20 |  3 |   4 |   8 |   1 |   2 |   7 |  18 |  11 |   6 |  13 |  14 |  12 |
| A5  | 12 | 17 | 19 | 15 |  7 | 16 |  2 | 20 |  3 |   9 |  13 |   1 |   4 |   5 |  18 |  11 |   6 |   8 |  10 |  14 |
| A6  | 10 | 15 | 19 | 16 |  8 | 18 |  6 | 20 |  3 |   7 |  11 |   1 |   2 |   4 |  17 |   9 |   5 |  13 |  12 |  14 |
| A7  |  9 | 16 | 19 | 17 | 10 | 15 |  5 | 20 |  3 |   8 |  11 |   1 |   2 |   6 |  18 |   7 |   4 |  14 |  12 |  13 |
| A8  | 14 | 18 | 20 | 19 | 11 | 15 |  6 | 17 |  4 |   3 |  10 |   1 |   2 |   7 |  16 |   8 |   5 |  12 |   9 |  13 |
| A9  |  8 | 16 | 18 | 19 | 12 | 13 |  6 | 20 |  5 |   3 |   7 |   1 |   4 |   2 |  17 |  10 |   9 |  15 |  14 |  11 |
| A10 |  7 | 17 | 19 | 18 |  9 | 15 |  5 | 20 |  3 |  10 |  11 |   1 |   2 |   6 |  16 |   8 |   4 |  13 |  12 |  14 |
| A11 | 12 | 16 | 19 | 15 | 13 | 18 |  7 | 20 |  3 |   5 |  11 |   1 |   2 |   6 |  17 |  10 |   4 |  14 |   8 |   9 |
| A12 | 14 | 15 | 19 | 16 | 12 | 18 |  8 | 20 |  3 |   4 |   9 |   1 |   2 |   7 |  17 |   6 |   5 |  13 |  10 |  11 |

Example dataset `potato_visual`.

## Algorithm Tuning

The `compute_mallows` function is the workhorse of `BayesMallows`. It
runs the Metropolis-Hastings algorithm and returns the posterior
distribution of the scale parameter \(\alpha\) and the latent ranks
\(\rho\) of the Mallows model. To see all its arguments, please run
`?compute_mallows` in the console.

We start by using all the default values of the parameters, so we only
need to supply the matrix of ranked items. We use the `potato_visual`
data printed above.

``` r
model_fit <- compute_mallows(potato_visual)
```

The argument returned is a list object of class `BayesMallows`, which
contains a whole lot of information about the MCMC run.

``` r
str(model_fit)
#> List of 16
#>  $ rho              : num [1:20, 1:3000] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:20] "P1" "P2" "P3" "P4" ...
#>   .. ..$ : NULL
#>  $ rho_acceptance   : num [1:3000, 1] 1 1 1 1 1 1 0 1 1 1 ...
#>  $ alpha            : num [1:3000, 1] 1 1 1 1 0.937 ...
#>  $ alpha_acceptance : num [1:3000, 1] 1 0 0 0 1 1 0 1 0 1 ...
#>  $ any_missing      : logi FALSE
#>  $ augpair          : logi FALSE
#>  $ metric           : chr "footrule"
#>  $ lambda           : num 0.1
#>  $ nmc              : int 3000
#>  $ n_items          : int 20
#>  $ n_assessors      : int 12
#>  $ alpha_jump       : int 1
#>  $ thinning         : int 1
#>  $ L                : int 4
#>  $ sd_alpha         : num 0.1
#>  $ aug_diag_thinning: int 100
#>  - attr(*, "class")= chr "BayesMallows"
```

The function `assess_convergence` produces plots for visual convergence
assessment. We start by studing \(\alpha\), which is the default. The
plot is shown below, and looks good enough, at least to begin with.

``` r
assess_convergence(model_fit)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-147-1.png)<!-- -->

Next, we study the convergence of \(\rho\). To avoid too complicated
plots, we pick 5 items to plot. Again, you can read more about this
function by running `assess_convergence` in the console.

``` r
assess_convergence(model_fit, type = "rho", items = 1:5)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-148-1.png)<!-- -->

When the name of the items have been specified in to column names of
`R`, we can also provide these to the `items` argument. The line below
plots potatoes 16 through
20.

``` r
assess_convergence(model_fit, type = "rho", items = c("P16", "P17", "P18", "P19", "P20"))
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-149-1.png)<!-- -->

Based on these plots, it looks like the algorithm starts to converge
after around 1000 iterations. Discarding the first 2000 iterations as
burn-in hence seems like a safe choice.

## Posterior Distributions

Once we are confident that the algorithm parameters are reasonable, we
can study the posterior distributions of the model parameters using the
generic function `plot.BayesMallows`. (Summary and print methods will be
added later.)

### Scale Parameter \(\alpha\)

With a burnin of 2000, the original `model_fit` object from the previous
subsection has 1000 MCMC samples. The default parameter of
`plot.BayesMallows` is \(alpha\), so we can study the posterior
distribution with the simple statement below.

``` r
plot(model_fit, burnin = 2000)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-150-1.png)<!-- -->

We see that the posterior distribution is a bit bumpy. Since this might
be due to randomness in the MCMC algorithm, we do another run with 10
times as many samples.

``` r
model_fit_big <- compute_mallows(potato_visual, nmc = 1e4 + 2000)
```

Next, we plot the posterior of \(\alpha\) from this run:

``` r
plot(model_fit_big, burnin = 2000)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-152-1.png)<!-- -->

This estimate of the posterior density of \(\alpha\) does not have the
bump to the right, so we clearly needed more samples. You can also try
again with, e.g., `1e5` samples, and see if the bimodality of this last
plot disappears as well (it does\!). We did not include that here,
because package vignettes should build rather fast.

The objects returned from `compute_mallows` contain the full MCMC
samples, so we remove `model_fit_big` before going on.

``` r
rm(model_fit_big)
```

### Latent Ranks \(\rho\)

Obtaining posterior samples from \(\rho\) is in general harder than for
\(\alpha\). Some items tend to be very sticky. We start by plotting the
`model_fit` object from above, with 3000 iterations, discarding the
first 2000 as burn-in. We now have to tell `plot.BayesMallows` that we
want a plot of `type = "rho"` and all the items. This gives us posterior
the posterior density of all the items.

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-154-1.png)<!-- -->

Some of the histograms in the plot above seem unreasonable peaked, so we
would like to try again with more samples.

#### Jumping of \(\alpha\)

Updating \(\alpha\) in each step may be time consuming, so we set
`alpha_jump = 10`. This implies that \(\alpha\) is updated and saved
only every 10th iteration on \(\rho\). Note that this is not thinning,
since \(\alpha\) is saved each time it is sampled. Before computing the
posterior distribution with this parameter, we do a new convergence
assessment, and hence generate a `BayesMallows` object with a small
number of iterations.

``` r
test_run <- compute_mallows(potato_visual, nmc = 10000, alpha_jump = 10)
```

The trace indicates convergence after around 200 iterations of
\(\alpha\), i.e., after 2000 Monte Carlo samples.

``` r
assess_convergence(test_run, type = "alpha")
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-156-1.png)<!-- -->

The convergence plot for \(\rho\) agrees that the MCMC algorithms seems
to have converged after 2000 iterations.

``` r
assess_convergence(test_run, type = "rho", items = 1:5)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-157-1.png)<!-- -->

We delete the `test_run` object and try with 25,000
iterations.

``` r
rm(test_run)
```

``` r
model_fit <- compute_mallows(potato_visual, nmc = 25000 + 2000, alpha_jump = 10)
```

The posterior density of \(\alpha\) looks similar to what it did above,
despite our use of `alpha_jump`.

``` r
plot(model_fit, burnin = 2000)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-160-1.png)<!-- -->

But the more interesting case is the latent ranks, which we plot below:

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-161-1.png)<!-- -->

In the plot of latent ranks above, potatoes 3, 6, 8, 9, 10, 12, 13, and
16 had no variation. In this plot, potatoes 3, 8, 12, and 13 have no
variation. Hence, adding more iterations seemed to help against
stickiness. In a real application, we would recommend running an even
larger sample.

#### Thinning

Saving a large number of iterations of \(\rho\) gets quite expensive, so
`compute_mallows` has a `thinning` parameter. It specifies that only
each `thinning`th iteration of \(\rho\) should be saved to memory. We
double the number of iterations, while setting `thinning = 2`. This
gives us the same number of posterior samples.

Please be careful with thinning. In this small data example it is
definitely wasteful\! Running the same number of iterations without
thinning always gives a better approximation of the posterior
distribution. Thinning might be useful when you need to run a large
number of iterations to explore the space of latent ranks, and the
latent ranks from all iterations do not fit in memory. (See, e.g.,
Gelman et al. (2004) for a discussion of thinning).

``` r
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 2)
```

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-163-1.png)<!-- -->

We can compare the computing time with or without thinning:

``` r
t1 <- Sys.time()
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 1)
t2 <- Sys.time()
(no_thinning <- t2 - t1)
#> Time difference of 0.5414681 secs
t3 <- Sys.time()
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 10)
t4 <- Sys.time()
(thinning <- t4 - t3)
#> Time difference of 0.321007 secs
```

With these data, and with 50,000 iterations, using thinning certainly
does not speed up the algorithm.

## Varying the Distance Metric

We can try to use the Kendall distance instead of the footrule distance.

``` r
model_fit <- compute_mallows(potato_visual, metric = "kendall",
                             nmc = 25000 + 2000, alpha_jump = 10)
```

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-166-1.png)<!-- -->

And we can use Spearman distance. In this case, since the number of
potatoes (20) is larger than the maximum number for which we have an
exact computation of the partition function (14), so precomputed
importance sampling estimates are used. This is handled automatically by
`compute_mallows`. Note that the posterior ranks are less peaked with
Spearman distance. This agrees with the results seen in Liu et al.
(2018).

``` r
model_fit <- compute_mallows(potato_visual, metric = "spearman",
                             nmc = 25000 + 2000, alpha_jump = 10)
```

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-168-1.png)<!-- -->

## Validation of Input

It is also worth pointing out that `compute_mallows` checks if the input
data are indeed ranks. Let us try this by manipulating the first row of
`potato_visual`, giving rank 1 to the first two items:

``` r
potato_modified <- potato_visual
potato_modified[1, 1:2] <- 1

model_fit <- compute_mallows(potato_modified)
#> Error in compute_mallows(potato_modified): Not valid permutation.
```

# Top-\(k\) Rankings

## Encoding of Missing Ranks

Now imagine that the assessors in the potato experiment were asked to
rank only the top-five heaviest potatoes. We generate these data by
retaining only ranks 5 or higher in `potato_visual`, setting the rest to
`NA`. (I use `dplyr::if_else` because it is safer than the `ifelse`
function of base R).

``` r
library(dplyr)
potato_top <- potato_visual * if_else(potato_visual > 5, NA_integer_, 1L)
```

Here is the resulting rank
matrix:

|     | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 | P10 | P11 | P12 | P13 | P14 | P15 | P16 | P17 | P18 | P19 | P20 |
| --- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| A1  | NA | NA | NA | NA | NA | NA |  4 | NA |  3 |   5 |  NA |   1 |   2 |  NA |  NA |  NA |  NA |  NA |  NA |  NA |
| A2  | NA | NA | NA | NA | NA | NA | NA | NA |  4 |   3 |  NA |   1 |   2 |  NA |  NA |  NA |   5 |  NA |  NA |  NA |
| A3  | NA | NA | NA | NA | NA | NA | NA | NA | NA |   3 |  NA |   2 |   1 |   4 |  NA |   5 |  NA |  NA |  NA |  NA |
| A4  | NA | NA | NA | NA | NA | NA |  5 | NA |  3 |   4 |  NA |   1 |   2 |  NA |  NA |  NA |  NA |  NA |  NA |  NA |
| A5  | NA | NA | NA | NA | NA | NA |  2 | NA |  3 |  NA |  NA |   1 |   4 |   5 |  NA |  NA |  NA |  NA |  NA |  NA |
| A6  | NA | NA | NA | NA | NA | NA | NA | NA |  3 |  NA |  NA |   1 |   2 |   4 |  NA |  NA |   5 |  NA |  NA |  NA |
| A7  | NA | NA | NA | NA | NA | NA |  5 | NA |  3 |  NA |  NA |   1 |   2 |  NA |  NA |  NA |   4 |  NA |  NA |  NA |
| A8  | NA | NA | NA | NA | NA | NA | NA | NA |  4 |   3 |  NA |   1 |   2 |  NA |  NA |  NA |   5 |  NA |  NA |  NA |
| A9  | NA | NA | NA | NA | NA | NA | NA | NA |  5 |   3 |  NA |   1 |   4 |   2 |  NA |  NA |  NA |  NA |  NA |  NA |
| A10 | NA | NA | NA | NA | NA | NA |  5 | NA |  3 |  NA |  NA |   1 |   2 |  NA |  NA |  NA |   4 |  NA |  NA |  NA |
| A11 | NA | NA | NA | NA | NA | NA | NA | NA |  3 |   5 |  NA |   1 |   2 |  NA |  NA |  NA |   4 |  NA |  NA |  NA |
| A12 | NA | NA | NA | NA | NA | NA | NA | NA |  3 |   4 |  NA |   1 |   2 |  NA |  NA |  NA |   5 |  NA |  NA |  NA |

Example dataset potato\_top.

In Vitelli et al. (2018) it is shown that the unranked items do not
effect the MAP estimates of the ranked items in this top-k setting. In
this case, there are 8 potatoes which have been ranked, and so the
unranked potatoes should have uniform posterior distributions between 9
and 20. However, arriving at these uniform posteriors require a large
number of MCMC iterations, so we instead remove these items:

``` r
item_ranked <- apply(potato_top, 2, function(x) !all(is.na(x)))
potato_top <- potato_top[, item_ranked, drop = FALSE]
```

We are now left with this 12 by 8 matrix:

|     | P7 | P9 | P10 | P12 | P13 | P14 | P16 | P17 |
| --- | -: | -: | --: | --: | --: | --: | --: | --: |
| A1  |  4 |  3 |   5 |   1 |   2 |  NA |  NA |  NA |
| A2  | NA |  4 |   3 |   1 |   2 |  NA |  NA |   5 |
| A3  | NA | NA |   3 |   2 |   1 |   4 |   5 |  NA |
| A4  |  5 |  3 |   4 |   1 |   2 |  NA |  NA |  NA |
| A5  |  2 |  3 |  NA |   1 |   4 |   5 |  NA |  NA |
| A6  | NA |  3 |  NA |   1 |   2 |   4 |  NA |   5 |
| A7  |  5 |  3 |  NA |   1 |   2 |  NA |  NA |   4 |
| A8  | NA |  4 |   3 |   1 |   2 |  NA |  NA |   5 |
| A9  | NA |  5 |   3 |   1 |   4 |   2 |  NA |  NA |
| A10 |  5 |  3 |  NA |   1 |   2 |  NA |  NA |   4 |
| A11 | NA |  3 |   5 |   1 |   2 |  NA |  NA |   4 |
| A12 | NA |  3 |   4 |   1 |   2 |  NA |  NA |   5 |

Example dataset `potato_top`.

## Metropolis-Hastings Algorithm with Missing Ranks

The `compute_mallows` function automatically recognizes the `NA` values
as missing ranks, and augments the data, as described in Section 4.1 of
Vitelli et al. (2018). Let us try:

``` r
model_fit <- compute_mallows(potato_top)
```

Looking at the returned object, we see that `any_missing` is `TRUE`, so
`compute_mallows` has correctly detected that there are missing values.

``` r
str(model_fit)
#> List of 17
#>  $ rho              : num [1:8, 1:3000] 1 2 3 4 5 6 7 8 1 2 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:8] "P7" "P9" "P10" "P12" ...
#>   .. ..$ : NULL
#>  $ rho_acceptance   : num [1:3000, 1] 1 0 0 1 1 1 0 0 0 1 ...
#>  $ alpha            : num [1:3000, 1] 1 1.1 1.35 1.38 1.37 ...
#>  $ alpha_acceptance : num [1:3000, 1] 1 1 1 1 1 1 1 1 0 1 ...
#>  $ any_missing      : logi TRUE
#>  $ augpair          : logi FALSE
#>  $ aug_acceptance   : num [1:12, 1:30] 0.6 0.76 0.99 0.71 0.8 0.82 0.64 0.67 0.75 0.63 ...
#>  $ metric           : chr "footrule"
#>  $ lambda           : num 0.1
#>  $ nmc              : int 3000
#>  $ n_items          : int 8
#>  $ n_assessors      : int 12
#>  $ alpha_jump       : int 1
#>  $ thinning         : int 1
#>  $ L                : int 1
#>  $ sd_alpha         : num 0.1
#>  $ aug_diag_thinning: int 100
#>  - attr(*, "class")= chr "BayesMallows"
```

## Algorithm Tuning

### Convergence of Augmentation

When data augmentation is used, `compute_mallows` saves the average
acceptance rates for each iteration interval of size
`aug_diag_thinning`, which is an optional argument default value 100. We
can specify to `assess_convergence` that we want to study the
convergence of the data augmentation.

``` r
assess_convergence(model_fit, type = "augmentation")
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-176-1.png)<!-- -->

If we want to take a closer look at some of the curves, we can specify
this with the `assessors` argument. Let’s look at assessor 3, which
seems to have a very high acceptance rate. First, we do two other runs
and compare them. These will differ due to the stochastic MCMC
algorithm. We also increase the number of iterations to 10,000.

``` r
model_fit1 <- compute_mallows(potato_top, nmc = 10000)
model_fit2 <- compute_mallows(potato_top, nmc = 10000)
```

The trace of acceptance for assessor three for the two independent runs
are shown in the plots below. Although high, these are absolutely
acceptable.

``` r
assess_convergence(model_fit1, type = "augmentation", assessors = 3)
assess_convergence(model_fit2, type = "augmentation", assessors = 3)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-178-1.png)![](inst/vignette_files/figure-gfm/unnamed-chunk-178-2.png)

We should also check that \(\alpha\) and \(\rho\) show good convergence
behavior, as before. We use the model object `model_fit1` for this, and
delete the other ones that we have created so far.

``` r
rm(model_fit, model_fit2)
```

### Convergence of \(\alpha\)

The convergence of \(\alpha\) looks very good.

``` r
assess_convergence(model_fit1, type = "alpha")
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-180-1.png)<!-- -->

### Convergence of \(\rho\)

The latent ranks also seem to converge.

``` r
assess_convergence(model_fit1, type = "rho", items = 1:8)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-181-1.png)<!-- -->

Based on the analysis in this section, we assume conservatively that the
MCMC algorithm reaches convergence after 10,000 iterations.

Before going on, we delete the last model object:

``` r
rm(model_fit1)
```

## Posterior Distributions

We now run `compute_mallows` bit longer, to obtain 20,000 samples from
the posterior distribution. There is no need for thinning in this case,
since the data fit well into memory.

``` r
model_fit <- compute_mallows(potato_top, nmc = 1e4 + 2e4)
```

Here is the posterior distribution of the scale parameter:

``` r
plot(model_fit, burnin = 1e4)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-184-1.png)<!-- -->

And the posterior distribution of the latent ranks:

``` r
plot(model_fit, burnin = 1e4, type = "rho", items = 1:8)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-185-1.png)<!-- -->

## Other Distance Measures

Like for the complete ranks, we can vary the distance measure used in
the Mallows model. We can try with a Spearman
model:

``` r
model_fit <- compute_mallows(potato_top, nmc = 1e4 + 2e4, metric = "spearman")
```

As for the full ranks described in the intro vignette, the posterior
uncertainty is higher with the Spearman distance.

``` r
plot(model_fit, burnin = 1e4, type = "rho", items = 1:8)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-187-1.png)<!-- -->

# Ranks Missing at Random

If the ranks are missing at random, we cannot remove the unranked items
as we did for top-\(k\) rankings above. Let us assume that 10 % of the
data in `potato_visual` have disappeared due to a disk failure. We
generate these in the code chunk below:

``` r
missing_indicator <- if_else(
  runif(nrow(potato_visual) * ncol(potato_visual)) < 0.1,
                            NA_real_, 1)
potato_missing <- potato_visual * missing_indicator
```

The data now look like the
following:

|     | P1 | P2 | P3 | P4 | P5 | P6 | P7 | P8 | P9 | P10 | P11 | P12 | P13 | P14 | P15 | P16 | P17 | P18 | P19 | P20 |
| --- | -: | -: | -: | -: | -: | -: | -: | -: | -: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| A1  | 10 | 18 | NA | 15 |  6 | 16 |  4 | 20 |  3 |   5 |  12 |   1 |   2 |  NA |  17 |  NA |   7 |  14 |  13 |  11 |
| A2  | 10 | 18 | 19 | 17 | 11 | 15 |  6 | 20 |  4 |   3 |  13 |   1 |   2 |   7 |  16 |   8 |   5 |  NA |   9 |  14 |
| A3  | 12 | 15 | 18 | 16 | NA | NA |  7 | NA |  6 |   3 |   8 |  NA |   1 |   4 |  19 |   5 |   9 |  14 |  10 |  17 |
| A4  |  9 | 17 | NA | 16 | 10 | 15 |  5 | 20 |  3 |  NA |   8 |   1 |   2 |   7 |  18 |  11 |   6 |  13 |  14 |  12 |
| A5  | 12 | 17 | 19 | 15 |  7 | 16 |  2 | 20 |  3 |   9 |  13 |   1 |   4 |   5 |  18 |  11 |   6 |   8 |  10 |  14 |
| A6  | 10 | 15 | 19 | 16 |  8 | 18 |  6 | 20 |  3 |   7 |  11 |   1 |   2 |   4 |  NA |   9 |   5 |  13 |  12 |  14 |
| A7  |  9 | 16 | 19 | 17 | 10 | 15 |  5 | 20 |  3 |   8 |  11 |   1 |   2 |   6 |  18 |   7 |   4 |  14 |  12 |  13 |
| A8  | 14 | 18 | 20 | 19 | 11 | 15 |  6 | 17 |  4 |   3 |  10 |   1 |  NA |   7 |  16 |  NA |  NA |  12 |   9 |  13 |
| A9  |  8 | 16 | NA | 19 | 12 | 13 |  6 | 20 |  5 |   3 |   7 |   1 |   4 |   2 |  17 |  10 |   9 |  15 |  NA |  11 |
| A10 |  7 | 17 | NA | 18 |  9 | 15 |  5 | 20 |  3 |  10 |  11 |   1 |   2 |   6 |  16 |   8 |   4 |  13 |  12 |  14 |
| A11 | 12 | 16 | 19 | 15 | 13 | 18 |  7 | NA |  3 |   5 |  11 |   1 |   2 |   6 |  17 |  10 |   4 |  14 |   8 |  NA |
| A12 | 14 | 15 | 19 | 16 | 12 | NA | NA | 20 |  3 |  NA |   9 |   1 |   2 |  NA |  17 |   6 |   5 |  13 |  10 |  11 |

Example dataset `potato_missing`.

## Algorithm Tuning

We supply `potato_missing` to `compute_mallows` as before:

``` r
model_fit <- compute_mallows(potato_missing, nmc = 1e4)
```

The convergence of \(\alpha\) and \(\rho\) seem fine:

``` r
assess_convergence(model_fit)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-191-1.png)<!-- -->

``` r
assess_convergence(model_fit, type = "rho", items = 1:6)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-192-1.png)<!-- -->

We also plot the convergence of the data augmentation.

``` r
assess_convergence(model_fit, type = "augmentation")
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-193-1.png)<!-- -->

Note in the plot above that the proposed augmentation is always accepted
for assessors 2, 5, 6, 7, and 10. We can compare this to the number of
missing ranks for these assessors:

``` r
apply(potato_missing, 1, function(x) sum(is.na(x)))
#>  A1  A2  A3  A4  A5  A6  A7  A8  A9 A10 A11 A12 
#>   3   1   4   2   0   1   0   3   2   1   2   4
```

We see that these assessors have either 0 or 1 missing ranks. In the
former case, there is no augmentation at all, but we use convention that
the assessor’s complete data are accepted. In the latter case, when only
one rank is missing, the proposal distribution has a support set of size
1, which has probability 1 of being accepted.

## Posterior Distributions

Again, we can fit a final model, and plot the posterior histogram of the
latent ranks.

``` r
model_fit <- compute_mallows(potato_visual, nmc = 1e4)
```

``` r
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-196-1.png)<!-- -->

# Pairwise Preferences

Handling of pairwise preferences in the Mallows rank model is described
in Section 4.2 of Vitelli et al. (2018).

## Introduction

Let us start by considering a toy example with two assessors and five
items. Assessor 1 has stated a set of preferences
\[ \mathcal{B}_{1} = \left\{A_{1} \prec A_{2}, A_{2} \prec A_{5}, A_{4} \prec A_{5} \right\} \]
and assessor 2 has the set of preferences
\[ \mathcal{B}_{2} = \left\{ A_{1} \prec A_{2}, A_{2} \prec A_{3}, A_{3} \prec A_{4} \right\}. \]

### Data Model

Each time an assessor is asked to compare two objects, a measurement is
made. Therefore, in order to keep the data *tidy* (Wickham (2014)), we
define a dataframe in which each row corresponds to a pairwise
comparison. The columns (variables) are *the assessor*, *the bottom
item*, and *the top item*.

In the code snippet below, we define such a dataframe for the toy
example presented above:

``` r
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4
)

knitr::kable(pair_comp)
```

| assessor | bottom\_item | top\_item |
| -------: | -----------: | --------: |
|        1 |            1 |         2 |
|        1 |            2 |         5 |
|        1 |            4 |         5 |
|        2 |            1 |         2 |
|        2 |            2 |         3 |
|        2 |            3 |         4 |

### Transitive Closure

Next, we need to find the transitive closure for the set of pairwise
comparisons given by each user. `BayesMallows` comes with a function
`generate_transitive_closure` to do just this.

``` r
pair_comp_tc <- generate_transitive_closure(pair_comp)
```

As we can see, `pair_comp_tc` has an additional row containing the
relation \(A_{4} \prec A_{5}\) for assessor 1. For assessor 2,
\[\text{tc}(\mathcal{B}_{2}) = \mathcal{B}_{2} \cup \left\{ A_{1} \prec A_{3}, A_{1} \prec A_{4}, A_{2} \prec A_{4}\right\},\]
so three new rows have been added.

``` r
knitr::kable(pair_comp_tc)
```

| assessor | bottom\_item | top\_item |
| -------: | -----------: | --------: |
|        1 |            1 |         2 |
|        1 |            1 |         5 |
|        1 |            2 |         5 |
|        1 |            4 |         5 |
|        2 |            1 |         2 |
|        2 |            1 |         3 |
|        2 |            2 |         3 |
|        2 |            1 |         4 |
|        2 |            2 |         4 |
|        2 |            3 |         4 |

The dataframe returned by `generate_transitive_closure` inherits from
`tibble`, but has subclass `BayesMallowsTC`. The `compute_mallows`
function uses this information to ensure that the object provided has
been through the `generate_transitive_closure` function. If it has not,
`compute_mallows` will do it for us, but this may lead to additional
computing time when running several diagnostic runs and trying out
different parameters, since the transitive closure will be recomputed
each time.

``` r
class(pair_comp_tc)
#> [1] "BayesMallowsTC" "tbl_df"         "tbl"            "data.frame"
```

### Initial Ranking

We can also generate an initial ranking, consistent with the pairwise
comparisons. Again, `compute_mallows` will do it for us, but we may save
time by computing it once and for all before we starting running the
algorithms.

``` r
initial_ranking <- generate_initial_ranking(pair_comp_tc)
```

``` r
knitr::kable(initial_ranking, row.names = TRUE)
```

|   | 1 | 2 | 3 | 4 | 5 |
| - | -: | -: | -: | -: | -: |
| 1 | 1 | 3 | 4 | 2 | 5 |
| 2 | 1 | 2 | 3 | 4 | 5 |

### Mallows Model

Having generated the transitive closure of each assessor’s pairwise
preferences and the initial ranking, we can go on and use these as
inputs to the Mallows model.

``` r
model_fit <- compute_mallows(R = initial_ranking, P = pair_comp_tc)
```

The model object has `augpair` equal to `TRUE`, and contains the
`aug_acceptance` statisics.

``` r
str(model_fit)
#> List of 17
#>  $ rho              : num [1:5, 1:3000] 1 2 3 4 5 2 1 3 4 5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : chr [1:5] "1" "2" "3" "4" ...
#>   .. ..$ : NULL
#>  $ rho_acceptance   : num [1:3000, 1] 1 1 1 0 1 1 0 1 1 1 ...
#>  $ alpha            : num [1:3000, 1] 1 0.975 1.109 1.167 1.065 ...
#>  $ alpha_acceptance : num [1:3000, 1] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ any_missing      : logi FALSE
#>  $ augpair          : logi TRUE
#>  $ aug_acceptance   : num [1:2, 1:30] 0.97 0.98 0.98 0.99 0.99 0.96 0.87 0.95 0.97 0.99 ...
#>  $ metric           : chr "footrule"
#>  $ lambda           : num 0.1
#>  $ nmc              : int 3000
#>  $ n_items          : int 5
#>  $ n_assessors      : int 2
#>  $ alpha_jump       : int 1
#>  $ thinning         : int 1
#>  $ L                : int 1
#>  $ sd_alpha         : num 0.1
#>  $ aug_diag_thinning: int 100
#>  - attr(*, "class")= chr "BayesMallows"
```

We can study the acceptance rate of the proposed augmented ranks.

``` r
assess_convergence(model_fit, type = "augmentation")
```

![](inst/vignette_files/figure-gfm/unnamed-chunk-205-1.png)<!-- -->

Rather than digging deeper into this toy example, we go on with a real
application.

## Beach Preferences

The beach preference dataset is described in Section 6.2 of Vitelli et
al. (2018), and is available in the dataframe `beach_preferences` in
`BayesMallows`. In short, \(60\) assessors were each asked to perform a
random set of pairwise comparisons between pictures of \(15\) beaches.
The first few rows are listed
below.

``` r
knitr::kable(head(beach_preferences, 6), caption = "Example dataset `beach_preferences`")
```

| assessor | bottom\_item | top\_item |
| -------: | -----------: | --------: |
|        1 |            2 |        15 |
|        1 |            5 |         3 |
|        1 |           13 |         3 |
|        1 |            4 |         7 |
|        1 |            5 |        15 |
|        1 |           12 |         6 |

Example dataset `beach_preferences`

### Transitive Closures

We start by generating the transitive closure of the preferences.

``` r
beach_tc <- generate_transitive_closure(beach_preferences)
```

We can compare the dataframes before and after. We see that the number
of rows has been approximately doubled, and that `beach_tc` has subclass
`BayesMallowsTC` has it should.

``` r
str(beach_preferences)
#> Classes 'tbl_df', 'tbl' and 'data.frame':    1442 obs. of  3 variables:
#>  $ assessor   : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ bottom_item: num  2 5 13 4 5 12 14 7 8 14 ...
#>  $ top_item   : num  15 3 3 7 15 6 6 9 10 9 ...
str(beach_tc)
#> Classes 'BayesMallowsTC', 'tbl_df', 'tbl' and 'data.frame':  2921 obs. of  3 variables:
#>  $ assessor   : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ bottom_item: int  4 5 7 4 5 7 8 9 13 14 ...
#>  $ top_item   : int  1 1 1 3 3 3 3 3 3 3 ...
```

### Initial Ranking

Next, we generate an initial ranking.

``` r
beach_init_rank <- generate_initial_ranking(beach_tc)
```

We can also take a look at the first 6 rows in
it.

``` r
knitr::kable(head(beach_init_rank, 6))
```

| 1 |  2 |  3 |  4 |  5 |  6 |  7 |  8 |  9 | 10 | 11 | 12 | 13 | 14 | 15 |
| -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: | -: |
| 2 |  4 |  5 |  8 | 14 |  7 | 15 | 13 |  1 |  9 | 12 | 10 |  3 |  6 | 11 |
| 2 |  3 |  6 |  7 |  8 | 12 | 13 | 14 | 10 |  4 |  1 | 11 | 15 |  5 |  9 |
| 4 |  5 |  8 | 11 | 12 | 13 |  2 | 10 |  3 |  7 | 14 | 15 |  1 |  6 |  9 |
| 2 |  4 | 12 | 14 |  5 |  8 |  7 |  6 | 13 | 15 |  1 | 10 | 11 |  3 |  9 |
| 2 |  7 |  8 | 12 | 13 | 15 |  1 | 14 |  3 |  5 |  4 |  9 | 10 | 11 |  6 |
| 5 | 12 |  2 |  4 |  8 | 14 |  1 | 13 |  7 | 10 |  3 | 11 |  6 | 15 |  9 |

### Algorithm Tuning

We can now check the convergence, using the same tools as before.

``` r
test_run <- compute_mallows(R = beach_init_rank, P = beach_tc, nmc = 2)
#> Error in run_mcmc(R = t(R), nmc = nmc, pairwise = P, constrained = constrained, : randi(): incorrect distribution parameters: a must be less than b
```

As you can tell, there is a bug. And that is where I will start
tomorrow. \# Clustering of Assessors

# References

<div id="refs" class="references">

<div id="ref-asfaw2016">

Asfaw, D., V. Vitelli, O. Sorensen, E. Arjas, and A. Frigessi. 2016.
“Time‐varying Rankings with the Bayesian Mallows Model.” *Stat* 6 (1):
14–30. <https://onlinelibrary.wiley.com/doi/abs/10.1002/sta4.132>.

</div>

<div id="ref-crispino2018">

Crispino, M., E. Arjas, V. Vitelli, N. Barrett, and A. Frigessi. 2018.
“A Bayesian Mallows approach to non-transitive pair comparison data:
how human are sounds?” *Accepted for Publication in Annals of Applied
Statistics*.
<https://www.e-publications.org/ims/submission/AOAS/user/submissionFile/34225?confirm=6b603456>.

</div>

<div id="ref-gelman2004">

Gelman, Andrew, John B. Carlin, Hal S. Stern, and Donald B. Rubin. 2004.
*Bayesian Data Analysis*. 2nd ed. Chapman; Hall/CRC.

</div>

<div id="ref-liu2018">

Liu, Q., M. Crispino, I. Scheel, V. Vitelli, and A. Frigessi. 2018.
“Model-based learning from preference data.” *Manuscript*.

</div>

<div id="ref-vitelli2018">

Vitelli, V., O. Sorensen, M. Crispino, E. Arjas, and A. Frigessi. 2018.
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research* 18 (1): 1–49.
<http://jmlr.org/papers/v18/15-481.html>.

</div>

<div id="ref-wickham2014">

Wickham, Hadley. 2014. “Tidy Data.” *Journal of Statistical Software,
Articles* 59 (10): 1–23. <https://doi.org/10.18637/jss.v059.i10>.

</div>

</div>
