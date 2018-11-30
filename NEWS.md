# BayesMallows 0.1.1.9008
* Fixed bug with Cayley distance. For this distance, the computational shortcut on p. 8 of Vitelli et al. (2018), JMLR, does not work. However, it was still used. Now, Cayley distance is always computed with complete rank vectors.
* Fixed bug in the default argument `leap_size` to `compute_mallows`. It used to be `floor(n_items / 5)`, which evaluates to zero when `n_items <= 4`. Updated it to `max(1L, floor(n_items / 5))`.
* Added Hamming distance (`metric = "haming"`) as an option to `compute_mallows` and `sample_mallows`.

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
