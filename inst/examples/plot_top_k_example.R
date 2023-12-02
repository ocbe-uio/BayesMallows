set.seed(1)
# We use the example dataset with beach preferences. Se the documentation to
# compute_mallows for how to assess the convergence of the algorithm
# We need to save the augmented data, so setting this option to TRUE
model_fit <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(
    nmc = 1000, burnin = 500, save_aug = TRUE))
# By default, the probability of being top-3 is plotted
# The default plot gives the probability for each assessor
plot_top_k(model_fit)
# We can also plot for each item
plot_top_k(model_fit, plot_level = "item")
# We can also plot the probability of being top-5, for each item
plot_top_k(model_fit, k = 5)
# We get the underlying numbers with predict_top_k
probs <- predict_top_k(model_fit)
# To find all items ranked top-3 by assessors 1-3 with probability more than 80 %,
# we do
subset(probs, assessor %in% 1:3 & prob > 0.8)

# We can also plot for clusters
model_fit <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  model_options = set_model_options(n_clusters = 3),
  compute_options = set_compute_options(
    nmc = 1000, burnin = 500, save_aug = TRUE)
  )
# In the item level plot, we now get one column per cluster
plot_top_k(model_fit, plot_level = "item")
