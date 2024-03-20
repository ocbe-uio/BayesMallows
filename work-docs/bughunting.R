set.seed(1)
user_ids <- 1:12
dat0 <- potato_visual
dat0[] <- ifelse(runif(length(dat0)) > .2, NA_real_, dat0)

mod0 <- compute_mallows(
  data = setup_rank_data(rankings = dat0),
  compute_options = set_compute_options(burnin = 1000)
)

dat1 <- potato_visual
dat1 <- ifelse(is.na(dat0) & runif(length(dat1)) > .5, NA_real_, dat1)

mod1 <- update_mallows(
  model = mod0,
  new_data = setup_rank_data(rankings = dat1, user_ids = user_ids),
  compute_options = set_compute_options(aug_method = "pseudo")
)
for(i in 1:12) {
  expect_equal(
    as.numeric(mod1$augmented_rankings[, i, 1][!is.na(dat1[i, ])]),
    as.numeric(dat1[i, !is.na(dat1[i, ])])
  )
}

mod_bmm1 <- compute_mallows(
  data = setup_rank_data(rankings = dat1),
  compute_options = set_compute_options(burnin = 500)
)

dat2 <- potato_visual
dat2 <- ifelse(is.na(dat1) & runif(length(dat2)) > .2, NA_real_, dat2)

mod2 <- update_mallows(
  model = mod1,
  new_data = setup_rank_data(rankings = dat2, user_ids = user_ids)
)

for(i in 1:12) {
  expect_equal(
    as.numeric(mod2$augmented_rankings[, i, 1][!is.na(dat2[i, ])]),
    as.numeric(dat2[i, !is.na(dat2[i, ])])
  )
}
