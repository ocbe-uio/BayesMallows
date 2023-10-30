# The example dataset beach_preferences contains pairwise preferences of beach.
# We must first generate the transitive closure
beach_tc <- generate_transitive_closure(beach_preferences)

# Next, we generate an initial ranking
beach_init <- generate_initial_ranking(beach_tc)

# Look at the first few rows:
head(beach_init)

# We can add more informative column names in order
# to get nicer posterior plots:
colnames(beach_init) <- paste("Beach", seq(from = 1, to = ncol(beach_init), by = 1))
head(beach_init)

# By default, the algorithm for generating the initial ranking is deterministic.
# It is possible to randomly permute the unranked items with the argument shuffle_unranked,
# as demonstrated below. This algorithm is computationally efficient, but defaults to FALSE
# for backward compatibility.
set.seed(2233)
beach_init <- generate_initial_ranking(beach_tc, shuffle_unranked = TRUE)
head(beach_init)

# It is also possible to pick a random sample among all topological sorts.
# This requires first enumerating all possible sorts, and might hence be computationally
# demanding. Here is an example, where we reduce the data considerable to speed up computation.
small_tc <- beach_tc[beach_tc$assessor %in% 1:6 &
                       beach_tc$bottom_item %in% 1:4 & beach_tc$top_item %in% 1:4, ]
set.seed(123)
init_small <- generate_initial_ranking(tc = small_tc, n_items = 4, random = TRUE)
# Look at the initial rankings generated
init_small

# For this small dataset, we can also study the effect of setting shuffle_unranked=TRUE
# in more detail. We consider assessors 1 and 2 only.
# First is the deterministic ordering. This one is equal for each run.
generate_initial_ranking(tc = small_tc[small_tc$assessor %in% c(1, 2), ],
                         n_items = 4, shuffle_unranked = FALSE, random = FALSE)
# Next we shuffle the unranked, setting the seed for reproducibility.
# For assessor 1, item 2 is unranked, and by rerunning the code multiple times,
# we see that element (1, 2) indeed changes ranking randomly.
# For assessor 2, item 3 is unranked, and we can also see that this item changes
# ranking randomly when rerunning the function multiple times.
# The ranked items also change their ranking from one random realiziation to another,
# but their relative ordering is constant.
set.seed(123)
generate_initial_ranking(tc = small_tc[small_tc$assessor %in% c(1, 2), ],
                         n_items = 4, shuffle_unranked = TRUE, random = FALSE)



# We now give beach_init and beach_tc to compute_mallows. We tell compute_mallows
# to save the augmented data, in order to study the convergence.
model_fit <- compute_mallows(
  rankings = beach_init, preferences = beach_tc,
  compute_options = set_compute_options(nmc = 1000, save_aug = TRUE)
)

# We can study the acceptance rate of the augmented rankings
assess_convergence(model_fit, parameter = "Rtilde")

# We can also study the posterior distribution of the consensus rank of each beach
model_fit$burnin <- 500
plot(model_fit, parameter = "rho", items = 1:15)

# The computations can also be done in parallel
library(parallel)
cl <- makeCluster(2)
beach_tc <- generate_transitive_closure(beach_preferences, cl = cl)
beach_init <- generate_initial_ranking(beach_tc, cl = cl)
stopCluster(cl)

