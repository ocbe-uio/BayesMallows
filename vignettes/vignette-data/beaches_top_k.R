library(BayesMallows)
library(ggplot2)

beach_tc <- generate_transitive_closure(beach_preferences)
beach_init_rank <- generate_initial_ranking(beach_tc)

model_fit <- compute_mallows(
  rankings = beach_init_rank,
  preferences = beach_tc,
  nmc = 6000,
  save_augmented_data = TRUE
)

ggsave(filename = "./vignettes/vignette-data/top_k_beaches.png",
       plot_top_k(model_fit, burnin = 2000))

