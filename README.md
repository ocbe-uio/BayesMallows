
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesMallows

This package implements the Bayesian Mallows Model described in Vitelli
et al. (2018). At the moment, the Mallows model for complete data is
implemented, and the user can choose between footrule, Spearman, Cayley,
Kendall, or Hamming distance.

Future releases will include:

  - Clustering users with similar preferences (Vitelli et al. (2018)).

  - Handling missing ranks by imputation (Vitelli et al. (2018)).

  - Time-varying ranks (Asfaw et al. (2016)).

  - Non-transitive pairwise comparisons (Crispino et al. (2017)).

All feedback and suggestions are very welcome. Feel free to create a
[Pull Request](https://github.uio.no/oyss/BayesMallows/pulls) or an
[Issue](https://github.uio.no/oyss/BayesMallows/issues).

Please also check out the [Suggested
Roadmap](https://github.uio.no/oyss/BayesMallows/wiki/Roadmap).

To install the current development version of the package, you should
clone or download this repository, and then open `BayesMallows.Rproj` in
RStudio and click **Build** and then **Install and Restart** on the top
menu. To get started using the package, take a look at the examples in
the `assess_convergence` and `compute_mallows` functions:

``` r
library(BayesMallows)
?assess_convergence
?compute_mallows
```

## Function Overview

This table describes the high level functions informally. It is aimed at
giving an overview of the package. Check the documentation for arguments
and return
values.

| Function Name            | Description                                                                                                                                                                                           |
| ------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `get_partition_function` | Returns the partition function for the Mallows model, given the necessary parameters, including distance measure. This function should have reasonable defaults for how the computation is performed. |
| `get_rank_distance`      | Returns the distance between two rank vectors, in the metric specified.                                                                                                                               |
| `compute_mallows`        | Runs the MCMC algorithm and computes the posterior distribution. The return values must include some convergence diagnostics. I suggest giving the return object S3 class `BayesMallows`.             |
| `assess_convergence`     | Function for investigating the convergence of the MCM algorithm, and trying different parameters of the proposal distributions.                                                                       |
| `plot.BayesMallows`      | S3 method for plotting the posterior distribution. Note that the user only needs to call `plot`.                                                                                                      |
| `summary.BayesMallows`   | S3 method for summarizing the posterior distribution.                                                                                                                                                 |

## References

<div id="refs" class="references">

<div id="ref-asfaw2016">

Asfaw, D., V. Vitelli, O. Sorensen, E. Arjas, and A. Frigessi. 2016.
“Time‐varying Rankings with the Bayesian Mallows Model.” *Stat* 6 (1):
14–30.

</div>

<div id="ref-crispino2017">

Crispino, M., E. Arjas, V. Vitelli, N. Barrett, and A. Frigessi. 2017.
“A Bayesian Mallows approach to non-transitive pair comparison data:
how human are sounds?” *ArXiv E-Prints*.

</div>

<div id="ref-vitelli2018">

Vitelli, V., O. Sorensen, M. Crispino, E. Arjas, and A. Frigessi. 2018.
“Probabilistic Preference Learning with the Mallows Rank Model.” *J
Mach. Learn. Res.* 18 (1): 1–49.

</div>

</div>
