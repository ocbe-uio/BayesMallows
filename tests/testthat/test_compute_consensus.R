context("Testing compute_consensus")

b <- compute_mallows(preferences = beach_preferences, nmc = 500)
b$burnin <- 200
cp <- compute_consensus(b)
map <- compute_consensus(b, type = "MAP")

expect_true(inherits(cp, "data.frame"))
expect_true(inherits(map, "data.frame"))

