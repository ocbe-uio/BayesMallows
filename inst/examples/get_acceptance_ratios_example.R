set.seed(1)
mod <- compute_mallows(
  data = setup_rank_data(potato_visual),
  compute_options = set_compute_options(burnin = 200)
)

get_acceptance_ratios(mod)
