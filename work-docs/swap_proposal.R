library(BayesMallows)
library(patchwork)
library(microbenchmark)

mod1 <- compute_mallows(
  data = setup_rank_data(potato_visual),
  compute_options = set_compute_options(nmc = 10000, burnin = 1000)
  )

mod2 <- compute_mallows(
  data = setup_rank_data(potato_visual),
  compute_options = set_compute_options(
    nmc = 10000, burnin = 1000, rho_proposal = "swap")
)

assess_convergence(mod1) + assess_convergence(mod2)

assess_convergence(mod1, parameter = "rho", items = 1:5) +
assess_convergence(mod2, parameter = "rho", items = 1:5)

plot(mod1) + plot(mod2)

microbenchmark(
  compute_mallows(
    data = setup_rank_data(potato_visual)
  ),
  compute_mallows(
    data = setup_rank_data(potato_visual),
    compute_options = set_compute_options(rho_proposal = "swap")
  )
)

mod1 <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 5000, burnin = 1000)
)
mod2 <- compute_mallows(
  data = setup_rank_data(preferences = beach_preferences),
  compute_options = set_compute_options(nmc = 5000, burnin = 1000, rho_proposal = "swap")
)
assess_convergence(mod1) + assess_convergence(mod2)
assess_convergence(mod1, parameter = "rho", items = 1:5) +
  assess_convergence(mod2, parameter = "rho", items = 1:5)

plot(mod1) + plot(mod2)

microbenchmark(
  compute_mallows(
    data = setup_rank_data(preferences = beach_preferences)
  ),
  compute_mallows(
    data = setup_rank_data(preferences = beach_preferences),
    compute_options = set_compute_options(rho_proposal = "swap")
  )
)
