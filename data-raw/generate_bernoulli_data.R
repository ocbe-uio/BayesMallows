set.seed(1)
theta <- .1
orderings <- as.data.frame(apply(potato_visual, 1, create_ordering))

prefs <- lapply(
  seq_along(orderings), function(i) {
    perfect_preferences <- matrix(rev(orderings[[i]]), ncol = 2, byrow = TRUE)
    noisy_preferences <- t(apply(perfect_preferences, 1, function(x) {
      if(runif(1) > theta) {
        x
      } else {
        rev(x)
      }
    }))
    dimnames(noisy_preferences) <- list(NULL, c("bottom_item", "top_item"))
    cbind(assessor = i, noisy_preferences)
  })
bernoulli_data <- do.call(rbind, prefs)

usethis::use_data(bernoulli_data, overwrite = TRUE)
