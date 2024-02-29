devtools::load_all()

dat <- subset(beach_preferences, assessor %in% 1:3)

mod <- compute_mallows(
  data = setup_rank_data(preferences = dat),
  compute_options = set_compute_options(nmc = 10000, burnin = 3000)
)


for(i in 4:20) {
  print(i)
  dat <- subset(beach_preferences, assessor == i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(preferences = dat)
  )
}
