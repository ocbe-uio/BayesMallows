library(BayesMallows)
library(dplyr)

pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4
)

pair_comp_tc <- generate_transitive_closure(pair_comp)
initial_ranking <- generate_initial_ranking(pair_comp_tc)

model_fit <- compute_mallows(rankings = initial_ranking,
                             preferences = pair_comp_tc, save_augmented_data = TRUE)

# In this plot, all augmented ranks get fixed after <500 iterations.
# I would expect that item 3 was drifting for assessor 1, and item 5 was drifting
# for assess 2.
assess_convergence(model_fit, type = "augmentation")
