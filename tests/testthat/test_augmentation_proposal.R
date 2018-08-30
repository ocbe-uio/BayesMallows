# context("Testing proposal distribution for data augmentation with preference data")
#
# library(BayesMallows)
# # This file is written to test that proposal of augmented
# # ranks for the case of pairwise preferences, described on
# # pp. 21-22 of Vitelli et al., JMLR (2018)
#
# # I test by creating an example dataset, and then walking through
# # step-by-step
#
# # The function BayesMallows:::find_pairwise_limits(u, ordering, current_ranking)
# # does the job of finding l_j and r_j, given an item u. Not that this is for a single
# # assessor.
#
# # One can look at: args(BayesMallows:::find_pairwise_limits)
# #
# # We set up an example dataset to test this. Here is a set of pairwise comparisons:
# pair_comp <- dplyr::tribble(
#   ~assessor, ~bottom_item, ~top_item,
#   1, 1, 2,
#   1, 2, 5,
#   1, 4, 5,
#   2, 1, 2,
#   2, 2, 3,
#   2, 3, 4,
#   3, 1, 4
# )
# n_items <- 5
#
# # We generate the transitive closure:
# pair_comp_tc <- generate_transitive_closure(pair_comp)
#
# # We then create a list of constraints per assessor
# constraints <- BayesMallows:::generate_constraints(pair_comp_tc, n_items)
#
# # The constraints list is used internally and compute_mallows(), and hence the generate_constraints() function
# # needs to be accessed with :::.
# # The contraints list contains one element per assessor.
# # For each assessor, i.e., list element, there are 3 elements. You can see them by running
# # For assessor 1:
# # R> j <- 1
# # R> names(constraints[[j]])
#
# # - constrained_items: a vector of the unique items that are constrained for assessor j
# # - items_above: a list. The ith element of items_above is the set of unique items that are preferred to item i
# # - items_below: a list. The ith element of items_below is the set of unique items that are disfavored to item i
#
# # You can run the following examples to understand:
# # First, find the constrained items
# # j <- 1
# # R> constraints[[j]][1]
# # Next, find the items above each item
# # R> constraints[[j]][2]
# # Finally, find the items below each item
# # R> constraints[[j]][3]
#
# # We set an intial ranking manually
# init_rank <- matrix(c(5, 5, 4, 2, 4, 1, 4, 3, 5, 3, 2, 2, 1, 1, 3), nrow = 3)
#
#
# # We are now ready to test the function for left and right limits
# # We cycle manually through u = 1,...,5 and j = 1, 2.
#
# #test_that("l_j and r_j of augmentation proposal are correct",{
#   #### Consider assessor 1 first
#   # Note that \mathcal{A}_{j} = \{1, 2, 4, 5\}
#   j <- 1
#
#   # Assume we draw the random number u = 1
#   # Now A_{j} \in \mathcal{A}_{j}.
#   u <- 1
#
#   # l_{j} is the maximum of the rankings of the elements that
#   # are preferred to item 1. According to constraints[[j]][["items_above"]][[u]],
#   # elements 2 and 5 are preferred to item 1. The maximum
#   # of their rankings is 2, so we expect l_{j} = 2.
#   # r_{j} is the minimum of the rankings of the elements that
#   # are disfavored to item 1. Since no items are disfavored to
#   # item 1, the set is empty, and we expect r_{j} = n_items + 1.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(3, n_items + 1)
#     )
#
#   # Next, assume we draw u = 2
#   u <- 2
#
#   # Item 5 is preferred to item 2. Item 5 has rank 1, so we expect l_{j} = 1.
#   # Items 1 and 4 are disfavored to item 2. The have rank 5 and 3, so we expect
#   # r_{j} = 3.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(1, 3)
#   )
#
#   # Next, assume we draw u = 3.
#   u <- 3
#
#   # Item 3 is not constrained, so we expect l_{j} = 0
#   # and r_{j} = n_items + 1.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, n_items + 1)
#   )
#
#   # Next, assume we draw u = 4
#   u <- 4
#
#   # Items 2 and 5 are preferred to 4, so we expect l_{j} = max(2, 1) = 2.
#   # Item 1 is disfavored to 4, so we expect r_{j} = min(5) = 5.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(2, 5)
#   )
#
#   # Next, assume we draw u = 5
#   u <- 5
#
#   # No items are preferred to item 5, so we assume l_{j} = 0.
#   # Items 1, 2, and 4 are disfavored to 5, so we assume r_{j} = min(5, 2, 3) = 2.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, 2)
#   )
#
#   #### Next, consider assessor 2
#   j <- 2
#
#   # Note that \mathcal{A}_{j} = \{1, 2, 3, 4\}
#
#   # Assume we draw the random number u = 1
#   # Now A_{j} \in \mathcal{A}_{j}.
#   u <- 1
#
#   # Items 2, 3, and 4 are preferred to 1, so we expect l_{j} = max(4, 3, 2) = 4.
#   # No items are disfavored to 1, so we expect r_{j} = n_items + 1 = 6
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(4, n_items + 1)
#   )
#
#   # Next, assume we draw u = 2
#   u <- 2
#
#   # Items 3 and 4 are preferred to 2, so we expect l_{j} = max(3, 2) = 3.
#   # Item 1 is disfavored to 2, so we expect r_{j} = min(5) = 5.
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(3, 5)
#   )
#
#   # Next, assume we draw u = 3.
#   u <- 3
#
#   # Item 4 is preferred to 3, so we expect l_{j} = max(2) = 2
#   # Items 1 and 2 are disfavored to 3, so we expect r_{j} = min(5, 4) = 4
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(2, 4)
#   )
#
#   # Next, assume we draw u = 4
#   u <- 4
#
#   # No items are preferred to 4, so we expect l_{j} = 0.
#   # Items 1, 2, and 3 are disfavored to 4, so we expect r_{j} = min(5, 4, 3) = 3
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, 3)
#   )
#
#   # Next, assume we draw u = 5
#   u <- 5
#
#   # Item 5 is not constrained, so we assume l_{j} = 0 and r_{j} = n_items + 1 = 6
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, n_items + 1)
#   )
#
#   #### Next, consider assessor 3
#   j <- 3
#
#   # Note that \mathcal{A}_{j} = \{1, 4\}
#
#   # Assume we draw the random number u = 1
#   # Now A_{j} \in \mathcal{A}_{j}.
#   u <- 1
#
#   # Item 4 is preferred to item 1, so we expect l_{j} = max(2) = 2
#   # No items are disfavored to item 1, so we expect r_{j} = n_items + 1 = 6
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(2, n_items + 1)
#   )
#
#   # Next, assume we draw u = 2
#   u <- 2
#
#   # Item 2 is not constrained, so we expect (0, n_items + 1)
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, n_items + 1)
#   )
#
#   # Next, assume we draw u = 3.
#   u <- 3
#
#   # Item 3 is not constrained, so we expect (0, n_items + 1)
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, n_items + 1)
#   )
#
#   # Next, assume we draw u = 4
#   u <- 4
#
#   # No items are preferred to 4, so we expect l_{j} = 0.
#   # Item 1 is disfavored to item 4, so we expect r_{j} = min(4) = 4
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, 4)
#   )
#
#   # Next, assume we draw u = 5
#   u <- 5
#
#   # Item 5 is not constrained, so we assume l_{j} = 0 and r_{j} = n_items + 1 = 6
#   expect_equal(
#     c(BayesMallows:::find_pairwise_limits(u, linear_ordering[[j]], init_rank[j, ])),
#     c(0, n_items + 1)
#   )
#
#           }
# )
#
# # Having passed the first test, we now go on to check that the proposed ranks
# # are consistent with the linear ordering. Under the hood, the find_pairwise_limits
# # function used in the last test is used to compute l_{j} and r_{j}
#
# # The function BayesMallows:::check_pairwise_augmentation takes an initial ranking
# # and a linear ordering, and proposes new rankings by using the full modified
# # leap-and-shift procedure. We hence use this to test that all proposals are consistent.
#
# # We test by always accepting the proposal, and running for a while, verifying
# # that inconsistent ranks are never proposed.
#
# # The function validate_rank_ordering is part of the BayesMallows package.
# proposal <- init_rank
# test_that("rank proposals are consistent with ordering", {
#
#   # Note that because the C++ code uses a transposed version of the rankings,
#   # we need to transpose both the argument and the result.
#   for(i in 1:100){
#     proposal <- t(BayesMallows:::check_pairwise_augmentation(t(proposal), linear_ordering))
#     expect_equal(validate_rank_ordering(proposal, linear_ordering), TRUE)
#   }
#
# })
#
# # Next, we run the same test on the beach data, which are more complex
#
# beach_tc <- generate_transitive_closure(beach_preferences)
#
# # We then create a linear ordering per assessor
# linear_ordering <- create_linear_ordering(beach_tc)
#
# # We generate a random initial ranking
# beach_init_rank <- generate_initial_ranking(beach_tc)
#
# # Set the proposal to the initial ranking
# proposal <- beach_init_rank
#
# test_that("beach dataset rank proposals are consistent with ordering", {
#
#   # Note that because the C++ code uses a transposed version of the rankings,
#   # we need to transpose both the argument and the result.
#   for(i in 1:100){
#     proposal <- t(BayesMallows:::check_pairwise_augmentation(t(proposal), linear_ordering))
#     expect_equal(validate_rank_ordering(proposal, linear_ordering), TRUE)
#   }
#
# })
