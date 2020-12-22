library(dplyr)
library(tidyr)
library(purrr)
# The first example uses full rankings in the potato_visual dataset, but we assume
# that each row in the data corresponds to between 100 and 500 assessors.
set.seed(1234)
# We start by generating random sample weights
weights <- sample(x = seq(from = 100L, to = 500L, by = 1L),
                  size = nrow(potato_visual), replace = TRUE)
# We also create a set of repeated indices, used to extend the matrix rows
repeated_indices <- unlist(map2(1:nrow(potato_visual), weights, ~ rep(.x, each = .y)))
# The potato_repeated matrix consists of all rows repeated corresponding to
# the number of assessors in the weights vector. This is how a large dataset
# would look like without using the weights argument
potato_repeated <- potato_visual[repeated_indices, ]

# We now first compute the Mallows model using weights
# This takes about 0.2 seconds
system.time({
  m_weights <- compute_mallows(rankings = potato_visual, weights = weights, nmc = 10000)
})
# Next we use the full ranking matrix
# This takes about 11.3 seconds, about 50 times longer!
\dontrun{
system.time({
  m_rep <- compute_mallows(rankings = potato_repeated, nmc = 10000)
})

  # We set the burnin to 2000 for both
  m_weights$burnin <- 2000
  m_rep$burnin <- 2000

  # Note that the MCMC algorithms did not run with the same
  # random number seeds in these two experiments, but still
  # the posterior distributions look similar
  plot(m_weights, burnin = 2000, "alpha")
  plot(m_rep, burnin = 2000, "alpha")

  plot(m_weights, burnin = 2000, "rho", items = 1:4)
  plot(m_rep, burnin = 2000, "rho", items = 1:4)
}

# Next we repeated the exercise with the pairwise preference data
# in the beach dataset. Note that we first must compute the
# transitive closure for each participant. If two participants
# have provided different preferences with identical transitive closure,
# then we can treat them as identical
beach_tc <- generate_transitive_closure(beach_preferences)
# Next, we confirm that each participant has a unique transitive closure
# We do this by sorting first by top_item and then by bottom_item,
# and then concatenating, whereupon we check how many participants there
# are for each unique concatenation
# This returns zero rows, so there are no participants with the same transitive closure
beach_tc %>%
  arrange(assessor, top_item, bottom_item) %>%
  group_by(assessor) %>%
  summarise(concat_ranks = paste(c(bottom_item, top_item), collapse = ","),
            .groups = "drop") %>%
  group_by(concat_ranks) %>%
  summarise(num_assessors = n_distinct(assessor), .groups = "drop") %>%
  filter(num_assessors > 1)

# We now illustrate the weighting procedure by assuming that there are
# more than one assessor per unique transitive closure. We generate a
# weights vector such that each unique transitive closure is repeated 1-4 times.
set.seed(9988)
weights <- sample(x = 1:4, size = length(unique(beach_preferences$assessor)), replace = TRUE)

# Next, we create a new hypthetical beach_preferences dataframe where each
# assessor is replicated 1-4 times
beach_pref_rep <- beach_preferences %>%
  mutate(new_assessor = map(weights[assessor], ~ 1:.x)) %>%
  unnest(cols = new_assessor) %>%
  mutate(assessor = paste(assessor, new_assessor, sep = ",")) %>%
  select(-new_assessor)

# We generate transitive closure for these preferences
beach_tc_rep <- generate_transitive_closure(beach_pref_rep)
# We can check that the number of unique assessors is now larger,
# and equal to the sum of weights
sum(weights)
length(unique(beach_tc_rep$assessor))

# We generate the initial rankings for the repeated and the "unrepeated"
# data
beach_rankings <- generate_initial_ranking(beach_tc, n_items = 15)
beach_rankings_rep <- generate_initial_ranking(beach_tc_rep, n_items = 15)

# We then run the Bayesian Mallows rank model, first for the
# unrepeated data with a weights argument. This takes about 1.9 seconds
system.time({
  model_fit_weights <- compute_mallows(rankings = beach_rankings,
                                       preferences = beach_tc,
                                       weights = weights,
                                       save_aug = TRUE,
                                       nmc = 10000)

})

# Next for the repeated data. This takes about 4.8 seconds.
system.time({
  model_fit_rep <- compute_mallows(rankings = beach_rankings_rep,
                                   preferences = beach_tc_rep,
                                   save_aug = TRUE,
                                   nmc = 10000)

})

# As demonstrated here, using a weights argument to exploit patterns in data
# where multiple assessors have given identical rankings or preferences, may
# lead to considerable speedup.
