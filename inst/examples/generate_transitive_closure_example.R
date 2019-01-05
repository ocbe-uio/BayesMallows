# Let us first consider a simple case with two assessors, where assessor 1
# prefers item 1 to item 2, and item 1 to item 5, while assessor 2 prefers
# item 3 to item 5. We then have the following dataframe of pairwise
# comparisons:
library(dplyr)
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 2, 1,
  1, 5, 1,
  2, 5, 3
)
# We then generate the transitive closure of these preferences:
(pair_comp_tc <- generate_transitive_closure(pair_comp))
# In this case, no additional relations we implied by the ones
# stated in pair_comp, so pair_comp_tc has exactly the same rows
# as pair_comp.

# Now assume that assessor 1 also preferred item 5 to item 3, and
# that assessor 2 preferred item 4 to item 3.
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 2, 1,
  1, 5, 1,
  1, 3, 5,
  2, 5, 3,
  2, 3, 4
)
# We generate the transitive closure again:
(pair_comp_tc <- generate_transitive_closure(pair_comp))
# We now have one implied relation for each assessor.
# For assessor 1, it is implied that 1 is preferred to 3.
# For assessor 2, it is implied that 4 is preferred to 5.

\dontrun{
  # If assessor 1 in addition preferred item 3 to item 1,
  # the preferences would not be consistent. This is not yet supported by compute_mallows,
  # so it causes an error message. It will be supported in a future release of the package.
  # First, we add the inconsistent row to pair_comp
  pair_comp <- bind_rows(pair_comp,
                         tibble(assessor = 1, bottom_item = 1, top_item = 3))

  # This causes an error message and prints out the problematic rankings:
  (pair_comp_tc <- generate_transitive_closure(pair_comp))
}


\dontrun{
  # The computations can also be done in parallel
  library(parallel)
  cl <- makeCluster(detectCores() - 1)
  beach_tc <- generate_transitive_closure(beach_preferences, cl = cl)
  stopCluster(cl)
}
