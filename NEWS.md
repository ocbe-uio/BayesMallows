# BayesMallows 0.2.0.9001
* Clusters are now `factor` variables sorted according to the cluster number. Hence, in plot legends, "Cluster 10" comes after "Cluster 9", rather than after "Cluster 1" which it used to do until now, because it was a `character`.
* `plot.BayesMallows` no longer contains print statements which forces display of plots. Instead plots are returned from the function. Using `p <- plot(fit)` hence does no longer display a plot, whereas using `plot(fit)` without assigning it to an object, displays a plot. Until now the plot was always shown for `rho` and `alpha`.

# BayesMallows 0.2.0.9000
* `compute_mallows` and `sample_mallows` now support Ulam distance, with argument `metric = "ulam"`.
* Slimmed down the vignette significantly, in order to avoid clang-UBSAN error caused by running the vignette (which was then again caused by `Rcpp`, cf. [this issue](https://github.com/RcppCore/Rcpp/issues/832)). The long vignette is no longer needed in any case, since all the functions are well documented with executable examples.

# BayesMallows 0.2.0
* New release on CRAN, which contains all the updates in 0.1.1, described below.

# BayesMallows 0.1.1.9009
* `Rankcluster` package has been removed from dependencies.

# BayesMallows 0.1.1.9008
* Fixed bug with Cayley distance. For this distance, the computational shortcut on p. 8 of Vitelli et al. (2018), JMLR, does not work. However, it was still used. Now, Cayley distance is always computed with complete rank vectors.
* Fixed bug in the default argument `leap_size` to `compute_mallows`. It used to be `floor(n_items / 5)`, which evaluates to zero when `n_items <= 4`. Updated it to `max(1L, floor(n_items / 5))`.
* Added Hamming distance (`metric = "hamming"`) as an option to `compute_mallows` and `sample_mallows`.

# BayesMallows 0.1.1.9007
* Updated `generate_initial_ranking`, `generate_transitive_closure`, and `sample_mallows` to avoid errors when package `tibble` version 2.0.0 is released. This update is purely internal.

# BayesMallows 0.1.1.9006
* Objects of class `BayesMallows` and `BayesMallowsMixtures` now have default print functions, hence avoiding excessive amounts of informations printed to the console if the user happens to write the name of such an object and press Return.
* `compute_mallows_mixtures` no longer sets `include_wcd = TRUE` by default. The user can choose this argument.
* `compute_mallows` has a new argument `save_clus`, which can be set to `FALSE` for not saving cluster assignments.

# BayesMallows 0.1.1.9005
* `assess_convergence` now automatically plots mixtures.
* `compute_mallows_mixtures` now returns an object of class `BayesMallowsMixtures`.

# BayesMallows 0.1.1.9004
* `assess_convergence` now adds prefix *Assessor* to plots when `parameter = "Rtilde"`.
* `predict_top_k` is now an exported function. Previously it was internal.

# BayesMallows 0.1.1.9003
* `compute_posterior_intervals` now has default `parameter = "alpha"`. Until now, this argument has had no default.
* Argument `type` to `plot.BayesMallows` and `assess_convergence` has been renamed to `parameter`, to be more consistent.

# BayesMallows 0.1.1.9002
* Argument `save_augment_data` to `compute_mallows` has been renamed to `save_aug`. 
* `compute_mallows` fills in implied ranks when an assessor has only one missing rank. This avoids unnecessary augmentation in MCMC.
* `generate_ranking` and `generate_ordering` now work with missing ranks.

# BayesMallows 0.1.1.9001
Argument `cluster_assignment_thinning` to `compute_mallows` has been renamed to `clus_thin`.

# BayesMallows 0.1.1.9000
Change the interface for computing consensus ranking. Now, CP and MAP consensus are both computed with the `compute_consensus` function, with argument `type` equal to either `"CP"` or `"MAP"`.
