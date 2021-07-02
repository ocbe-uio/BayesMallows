context("SMC pseudolikelihood functions")

################################################################################
# test for get_sample_probabilities
################################################################################
set.seed(101)
rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6

item_ordering = c(3,6,4,5)
partial_ranking = c(1,2,NA,NA,NA,NA)
remaining_set = c(3,4,5,6)
test_1 = get_sample_probabilities(rho_item_rank = rho[3], alpha, remaining_set_ranks = remaining_set, metric, n_items)
test_2 = get_sample_probabilities(rho_item_rank = rho[4], alpha, remaining_set_ranks = remaining_set, metric, n_items)

test_that('get_sample_probabilities outputs as expected', {
	expect_equivalent(test_1, c(0.3849370, 0.2758194, 0.1976332, 0.1416104), tol = 1e-6)
	expect_equivalent(test_2, c(0.2431822, 0.3393880, 0.2431822, 0.1742476), tol = 1e-6)
})

################################################################################
# test for calculate_forwards_probability and calculate_bacwards_probability
################################################################################
set.seed(101)
rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6

item_ordering = c(3,6,4,5)
partial_ranking = c(1,2,NA,NA,NA,NA)
remaining_set = c(3,4,5,6)

test_1_forward = calculate_forward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                               remaining_set = remaining_set, rho = rho, alpha = alpha,
                                               n_items = n_items, metric = metric)


current_ranking = c(1,2,6,5,4,3)
test_1_backward_a= calculate_backward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                                  current_ranking = current_ranking, remaining_set = remaining_set, rho = rho,
                                                  alpha = alpha, n_items = n_items, metric = metric)


current_ranking = test_1_forward$aug_ranking
test_1_backward_b = calculate_backward_probability(item_ordering = item_ordering, partial_ranking = partial_ranking,
                                                   current_ranking = current_ranking, remaining_set = remaining_set, rho = rho,
                                                   alpha = alpha, n_items = n_items, metric = metric)

test_that('calculations of forward and backward probabilities', {
    # print(test_1_forward)
    #$aug_ranking
    #[1] 1 2 3 5 4 6
    #$forward_prob
    #[1] 0.07205734

    # print(test_1_backward_a)
    # 0.01360987

    # print(test_1_backward_b)
    # 0.07205734

    # test_1_forward$forward_prob == test_1_backward_b
    # TRUE
})






############################################
# tests for M-H_aug_ranking_pseudo
###########################################



set.seed(101)
rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,4,5,6)

# test_1 = metropolis_hastings_aug_ranking_pseudo(alpha, rho, n_items, partial_ranking = R_obs, current_ranking = R_curr, metric) # FIXME: triggering error on get_rank_distance (#98)

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,NA,NA,NA)
# test_2 = metropolis_hastings_aug_ranking_pseudo(alpha, rho, n_items, partial_ranking = R_obs, current_ranking = R_curr, metric) # FIXME: triggering error on get_rank_distance (#98)


R_curr = c(1,2,6,5,4,3)
R_obs = c(1,2,NA,NA,NA,NA)
# test_3 =  metropolis_hastings_aug_ranking_pseudo(alpha, rho, n_items, partial_ranking = R_obs, current_ranking = R_curr, metric) # FIXME: triggering error on get_rank_distance (#98)

test_that('M-H aug ranking pseudo works', {
    # print(test_1)
    #$ranking
    #[1] 1 2 3 4 5 6

    # all(test_1 == R_curr) == TRUE
    # TRUE

    # print(test_2)
    # #$ranking
    # #[1] 1 2 3 5 4 6

    # all(test_2 == R_curr) == TRUE
    # #  should be FALSE

    # print(test_3)
    # #$ranking
    # #[1] 1 2 5 6 4 3

    # all(test_3 == R_curr) == TRUE
    # #  should be FALSE
})



#########################################################
### tests relating to the correction_kernel function
#########################################################


set.seed(101)
rho = c(1,2,3,4,5,6)
alpha = 2
metric = "footrule"
n_items= 6

R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,NA,NA,NA)

test_1 = correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)


R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,5,NA,NA)
test_2 = correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)


R_curr = c(1,2,3,4,5,6)
R_obs = c(1,2,3,4,5,6)
test_3 = correction_kernel_pseudo(R_curr, R_obs, rho, alpha, n_items, metric)

test_that('correction_kernel works', {
    # print(test_1)
    # #$ranking
    # #[1] 1 2 3 4 5 6
    # #
    # #$correction_prob
    # #[1] 1

    # all(test_1$ranking == R_curr) == TRUE
    # TRUE

    # print(test_2)
    # #$ranking
    # #[1] 1 2 3 5 4 6
    # #
    # #$correction_prob
    # #[1] 0.5


    # all(test_2$ranking == R_curr) == TRUE
    # #  should be FALSE

    # print(test_3)
    # #$ranking
    # #[1] 1 2 3 4 5 6
    # #
    # #$correction_prob
    # #[1] 1

    # all(test_3$ranking == R_curr) == TRUE
    # #  should be TRUE
})