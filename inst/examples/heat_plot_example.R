# Setting the number of Monte Carlo samples very low for the example to run fast.
# A real application should run much longer, and have a large burnin.
model_fit <- compute_mallows(
  potato_visual,
  compute_options = set_compute_options(nmc = 500),
  seed = 1)
model_fit$burnin <- 100

heat_plot(model_fit)

# Items are ordered along the horizontal axis according to the ordering
# returned by compute_consensus, whose default argument is type="CP".

heat_plot(model_fit, type = "MAP")
