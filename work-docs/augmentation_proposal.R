set.seed(1)
forward <- FALSE

rankings <- sample(10)
rho <- 1:10
unranked_items <- c(2, 4, 5, 8, 9, 10)
n_items <- length(rankings)
alpha <- 2
prob <- 1
while(length(unranked_items) > 0) {
  available_rankings <- rankings[unranked_items]
  item_to_rank <- unranked_items[[1]]
  log_numerator <- numeric()
  for(ll in 1:length(available_rankings)) {
    log_numerator[[ll]] <- -alpha / n_items * abs(rho[[item_to_rank]] - available_rankings[[ll]])
  }
  sample_probs <- exp(log_numerator) / sum(exp(log_numerator))

  if(forward) {
    ans <- rmultinom(1, seq_along(sample_probs), sample_probs)
    rankings[[item_to_rank]] <- available_rankings[ans == 1]
  }
  ranking_chosen <- which(rankings[[item_to_rank]] == available_rankings)
  prob <- prob * sample_probs[[ranking_chosen]]
  if(length(available_rankings) <= 1) break
  unranked_items <- unranked_items[-1]
  rankings[unranked_items] <-
    setdiff(available_rankings, available_rankings[[ranking_chosen]])
}

prob
