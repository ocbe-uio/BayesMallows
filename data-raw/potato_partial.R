library(BayesMallows)
set.seed(9999)
example_dataset <- potato_visual
n_users <- 12
n_items <- 20
potato_partial <- array(0, c(n_users, n_items, (n_items / 2 + 1)))
potato_partial[, , (n_items / 2 + 1)] <- potato_visual
tt <- 0
for (ii in (n_items - 1):(n_items / 2)) {
  tt <- tt + 1

  # set n_users line with one more NA
  example_dataset[example_dataset > ii] <- NA

  # set as new time stamp
  potato_partial[, , ((n_items / 2 + 1) - tt)] <- example_dataset
}

usethis::use_data(potato_partial, overwrite = TRUE)
