
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesMallows

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/BayesMallows)](https://cran.r-project.org/package=BayesMallows)
[![Build
Status](https://travis-ci.org/ocbe-uio/BayesMallows.svg?branch=master)](https://travis-ci.org/ocbe-uio/BayesMallows)
[![codecov](https://codecov.io/gh/ocbe-uio/BayesMallows/branch/master/graph/badge.svg)](https://codecov.io/gh/ocbe-uio/BayesMallows)

This package implements the Bayesian Mallows Model described in Vitelli
et al. (2018). The user can choose between footrule, Spearman, Cayley,
Hamming, Kendall, or Ulam distance.

The following features are currently implemented:

  - Complete data (Vitelli et al. (2018)).

  - Clustering users with similar preferences (Vitelli et al. (2018)).

  - Handling missing ranks by imputation (Vitelli et al. (2018)).

  - Handling transitive pairwise preferences by imputation (Vitelli et
    al. (2018)).

  - Estimating the partition function of the Mallows model using
    importance sampling (Vitelli et al. (2018)) or an asymptotic
    approximation (Mukherjee (2016)).

  - Non-transitive pairwise comparisons (Crispino et al. (2019)).

This includes any combination thereof, e.g., clustering assessors based
on pairwise preferences.

Future releases will include:

  - Time-varying ranks (Asfaw et al. (2016)).

  - Parallelization of Markov Chains.

All feedback and suggestions are very welcome.

## Installation

To install the current release, use

``` r
install.packages("BayesMallows")
```

To install the current development version, use

``` r
#install.packages("devtools")
devtools::install_github("ocbe-uio/BayesMallows")
```

## References

<div id="refs" class="references">

<div id="ref-asfaw2016">

Asfaw, Derbachew, Valeria Vitelli, Øystein Sørensen, Elja Arjas, and
Arnoldo Frigessi. 2016. “Time-Varying Rankings with the Bayesian Mallows
Model.” *Stat* 6 (1): 14–30. <https://doi.org/10.1002/sta4.132>.

</div>

<div id="ref-crispino2019">

Crispino, Marta, Elja Arjas, Valeria Vitelli, Natasha Barrett, and
Arnoldo Frigessi. 2019. “A Bayesian Mallows Approach to Nontransitive
Pair Comparison Data: How Human Are Sounds?” *The Annals of Applied
Statistics* 13 (1): 492–519. <https://doi.org/10.1214/18-aoas1203>.

</div>

<div id="ref-mukherjee2016">

Mukherjee, Sumit. 2016. “Estimation in Exponential Families on
Permutations.” *The Annals of Statistics* 44 (2): 853–75.
<https://doi.org/10.1214/15-aos1389>.

</div>

<div id="ref-vitelli2018">

Vitelli, V., O. Sorensen, M. Crispino, E. Arjas, and A. Frigessi. 2018.
“Probabilistic Preference Learning with the Mallows Rank Model.”
*Journal of Machine Learning Research* 18 (1): 1–49.
<http://jmlr.org/papers/v18/15-481.html>.

</div>

</div>
