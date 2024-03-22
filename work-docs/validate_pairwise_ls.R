# Validate pairwise leap-and-shift proposal
devtools::load_all()
library(tidyverse)
library(furrr)

dat <- subset(beach_preferences, assessor == 3 & bottom_item < 7 & top_item < 12)
n_items <- 15
dd <- setup_rank_data(preferences = dat, n_items = n_items)
constraints <- dd$constraints[[1]]

current_rank <- dd$rankings

plan("multisession")
rankings <- future_map_dfr(1:1000000, function(i) {
  u <- sample(n_items, 1)
  if(u %in% constraints$constrained_items) {
    ia <- constraints$items_above[[u]]
    ib <- constraints$items_below[[u]]
    if(length(ia) > 0) {
      l <- max(current_rank[ia]) + 1
    } else {
      l <- 1
    }
    if(length(ib) > 0) {
      r <- min(current_rank[ib]) - 1
    } else {
      r <- n_items
    }
  } else {
    l <- 1
    r <- n_items
  }
  support <- seq(from = l, to = r, by = 1)
  r_prime <- current_rank
  r_prime[[u]] <- sample(support, 1)

  for(j in seq_len(n_items)) {
    if(current_rank[[u]] < current_rank[[j]] && current_rank[[j]] <= r_prime[[u]]) {
      r_prime[[j]] <- current_rank[[j]] - 1
    } else if(current_rank[[u]] > current_rank[[j]] && current_rank[[j]] >= r_prime[[u]]) {
      r_prime[[j]] <- current_rank[[j]] + 1
    }
  }
  stopifnot(BayesMallows:::validate_permutation(r_prime))
  tibble(
    iteration = i,
    probability = 1 / length(support) / n_items,
    ranking = list(as.numeric(r_prime))
  )
}, .options = furrr_options(seed = 1L))

empirical_probs <- rankings %>%
  mutate(ranking = map_chr(ranking, paste, collapse = "")) %>%
  group_by(ranking) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(proportion = n / sum(n))

calculated_probs <- rankings %>%
  mutate(ranking = map_chr(ranking, paste, collapse = "")) %>%
  distinct(ranking, probability) %>%
  group_by(ranking) %>%
  summarise(
    probability = sum(probability), .groups = "drop"
  )

empirical_probs %>%
  full_join(calculated_probs, by = "ranking") %>%
  ggplot(aes(x = proportion, y = probability)) +
  geom_point() +
  geom_abline()
