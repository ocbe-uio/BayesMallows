# Package index

## Model estimation

Functions for estimating the Bayesian Mallows model.

- [`` `burnin<-`() ``](burnin-set.md) : Set the burnin
- [`burnin()`](burnin.md) : See the burnin
- [`compute_mallows()`](compute_mallows.md) : Preference Learning with
  the Mallows Rank Model
- [`compute_mallows_mixtures()`](compute_mallows_mixtures.md) : Compute
  Mixtures of Mallows Models
- [`compute_mallows_sequentially()`](compute_mallows_sequentially.md) :
  Estimate the Bayesian Mallows Model Sequentially
- [`sample_prior()`](sample_prior.md) : Sample from prior distribution
- [`update_mallows()`](update_mallows.md) : Update a Bayesian Mallows
  model with new users
- [`assess_convergence()`](assess_convergence.md) : Trace Plots from
  Metropolis-Hastings Algorithm

## Preparing model estimation

Functions to run prior to fitting a model.

- [`get_transitive_closure()`](get_transitive_closure.md) : Get
  transitive closure
- [`set_compute_options()`](set_compute_options.md) : Specify options
  for computation
- [`set_initial_values()`](set_initial_values.md) : Set initial values
  of scale parameter and modal ranking
- [`set_model_options()`](set_model_options.md) : Set options for
  Bayesian Mallows model
- [`set_priors()`](set_priors.md) : Set prior parameters for Bayesian
  Mallows model
- [`set_progress_report()`](set_progress_report.md) : Set progress
  report options for MCMC algorithm
- [`set_smc_options()`](set_smc_options.md) : Set SMC compute options
- [`setup_rank_data()`](setup_rank_data.md) : Setup rank data

## Posterior quantities

Functions for studying posterior distributions of model parameters.

- [`assign_cluster()`](assign_cluster.md) : Assign Assessors to Clusters
- [`compute_consensus()`](compute_consensus.md) : Compute Consensus
  Ranking
- [`compute_posterior_intervals()`](compute_posterior_intervals.md) :
  Compute Posterior Intervals
- [`get_acceptance_ratios()`](get_acceptance_ratios.md) : Get Acceptance
  Ratios
- [`heat_plot()`](heat_plot.md) : Heat plot of posterior probabilities
- [`plot(`*`<BayesMallows>`*`)`](plot.BayesMallows.md) : Plot Posterior
  Distributions
- [`plot(`*`<SMCMallows>`*`)`](plot.SMCMallows.md) : Plot SMC Posterior
  Distributions
- [`plot_elbow()`](plot_elbow.md) : Plot Within-Cluster Sum of Distances
- [`plot_top_k()`](plot_top_k.md) : Plot Top-k Rankings with Pairwise
  Preferences
- [`predict_top_k()`](predict_top_k.md) : Predict Top-k Rankings with
  Pairwise Preferences
- [`print(`*`<BayesMallows>`*`)`](print.BayesMallows.md)
  [`print(`*`<BayesMallowsMixtures>`*`)`](print.BayesMallows.md)
  [`print(`*`<SMCMallows>`*`)`](print.BayesMallows.md) : Print Method
  for BayesMallows Objects

## Rank functions

Various functions for sampling ranks and working with ranks.

- [`compute_expected_distance()`](compute_expected_distance.md) :
  Expected value of metrics under a Mallows rank model
- [`compute_observation_frequency()`](compute_observation_frequency.md)
  : Frequency distribution of the ranking sequences
- [`compute_rank_distance()`](compute_rank_distance.md) : Distance
  between a set of rankings and a given rank sequence
- [`create_ranking()`](create_ranking.md)
  [`create_ordering()`](create_ranking.md) : Convert between ranking and
  ordering.
- [`get_mallows_loglik()`](get_mallows_loglik.md) : Likelihood and
  log-likelihood evaluation for a Mallows mixture model
- [`sample_mallows()`](sample_mallows.md) : Random Samples from the
  Mallows Rank Model

## Partition functions

Tools related to computing or estimating the partition function of the
Mallows model with various distances.

- [`compute_exact_partition_function()`](compute_exact_partition_function.md)
  : Compute exact partition function
- [`estimate_partition_function()`](estimate_partition_function.md) :
  Estimate Partition Function
- [`get_cardinalities()`](get_cardinalities.md) : Get cardinalities for
  each distance

## Datasets

Example datasets included in the package.

- [`beach_preferences`](beach_preferences.md) : Beach preferences
- [`bernoulli_data`](bernoulli_data.md) : Simulated intransitive
  pairwise preferences
- [`cluster_data`](cluster_data.md) : Simulated clustering data
- [`potato_true_ranking`](potato_true_ranking.md) : True ranking of the
  weights of 20 potatoes.
- [`potato_visual`](potato_visual.md) : Potato weights assessed visually
- [`potato_weighing`](potato_weighing.md) : Potato weights assessed by
  hand
- [`sounds`](sounds.md) : Sounds data
- [`sushi_rankings`](sushi_rankings.md) : Sushi rankings
