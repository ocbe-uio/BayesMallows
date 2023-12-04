devtools::load_all()
set.seed(1)
theta <- 0.1
rankings <- potato_visual
dimnames(rankings) <- NULL

bernoulli_data <- expand.grid(
  comparison = 1:100,
  assessor = seq_len(nrow(rankings)),
  bottom_item = numeric(1),
  top_item = numeric(1)
  )

for(i in seq_len(nrow(bernoulli_data))) {
  r <- rankings[bernoulli_data$assessor[[i]], ]
  items_to_compare <- sample(seq_along(r), 2)

  a <- if(runif(1) > theta) {
      order(r[items_to_compare], decreasing = TRUE)
    } else {
      order(r[items_to_compare], decreasing = FALSE)
    }
  bernoulli_data[i, c("bottom_item", "top_item")] <- items_to_compare[a]
}

bernoulli_data$comparison <- NULL
usethis::use_data(bernoulli_data, overwrite = TRUE)
