library(BayesMallows)
library(tidyverse)
library(patchwork)

find_data <- function(k) {
  dat <- sushi_rankings[1:300, ]
  dat[dat > k] <- NA
  dat
}

mod <- compute_mallows(
  setup_rank_data(find_data(3)),
  compute_options = set_compute_options(burnin = 500))

plots <- list(plot(mod) + ggtitle(paste("k =", 3)))


for(i in 4:10) {
  print(i)
  mod <- update_mallows(
    model = mod,
    new_data = setup_rank_data(find_data(i))
  )
  plots <- c(plots, list(plot(mod) + ggtitle(paste("k = ", i))))
}

wrap_plots(plots)
