# Let us first consider a simple case with two assessors, where assessor 1
# prefers item 1 to item 2, and item 1 to item 5, while assessor 2 prefers
# item 3 to item 5. We then have the following dataframe of pairwise
# comparisons:
pair_comp <- data.frame(
  assessor = c(1, 1, 2),
  bottom_item = c(2, 5, 5),
  top_item = c(1, 1, 3))
# We then generate the transitive closure of these preferences:
(pair_comp_tc <- generate_transitive_closure(pair_comp))
# In this case, no additional relations we implied by the ones
# stated in pair_comp, so pair_comp_tc has exactly the same rows
# as pair_comp.

# Now assume that assessor 1 also preferred item 5 to item 3, and
# that assessor 2 preferred item 4 to item 3.
pair_comp <- data.frame(
  assessor = c(1, 1, 1, 2, 2),
  bottom_item = c(2, 5, 3, 5, 3),
  top_item = c(1, 1, 5, 3, 4))

# We generate the transitive closure again:
(pair_comp_tc <- generate_transitive_closure(pair_comp))
# We now have one implied relation for each assessor.
# For assessor 1, it is implied that 1 is preferred to 3.
# For assessor 2, it is implied that 4 is preferred to 5.

# The computations can also be done in parallel
library(parallel)
cl <- makeCluster(2)
beach_tc <- generate_transitive_closure(beach_preferences, cl = cl)
stopCluster(cl)

