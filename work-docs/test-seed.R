library(BayesMallows)
set.seed(1)
mod <- compute_mallows(setup_rank_data(potato_visual))
tail(mod$alpha$value)
