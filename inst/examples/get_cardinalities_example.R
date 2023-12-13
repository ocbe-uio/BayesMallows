# Extract the cardinalities for four items with footrule distance
n_items <- 4
dat <- get_cardinalities(n_items)
# Compute the partition function at alpha = 2
alpha <- 2
sum(dat$value * exp(-alpha / n_items * dat$distance))
#'
# We can confirm that it is correct by enumerating all possible combinations
all <- expand.grid(1:4, 1:4, 1:4, 1:4)
perms <- all[apply(all, 1, function(x) {length(unique(x)) == 4}),]
sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum(abs(x - 1:4)))))

# We do the same for the Spearman distance
dat <- get_cardinalities(n_items, metric = "spearman")
sum(dat$value * exp(-alpha / n_items * dat$distance))
#'
# We can confirm that it is correct by enumerating all possible combinations
sum(apply(perms, 1, function(x) exp(-alpha / n_items * sum((x - 1:4)^2))))
