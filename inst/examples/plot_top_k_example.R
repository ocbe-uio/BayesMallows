# We use the example dataset with beach preferences. Se the documentation to
# compute_mallows for how to assess the convergence of the algorithm
# We need to save the augmented data, so setting this option to TRUE
model_fit <- compute_mallows(preferences = beach_preferences,
                             save_augmented_data = TRUE)
# We set burnin = 1000
burnin <- 1000
# By default, the probability of being top-3 is plotted
plot_top_k(model_fit, burnin = burnin)
# We can also plot the probability of being top-5, for each item
plot_top_k(model_fit, burnin = burnin, k = 5)



