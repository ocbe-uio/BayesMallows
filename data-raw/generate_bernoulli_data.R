devtools::load_all()
set.seed(1)
theta <- 0.1
rankings <- potato_visual
dimnames(rankings) <- NULL

items <- expand.grid(1:20, 1:20)
items <- items[items[, 1] != items[, 2], ]
items <- items[items[, 1] < items[, 2], ]

bernoulli_data <- expand.grid(
  comp = seq_len(nrow(items)),
  assessor = seq_len(nrow(rankings)),
  bottom_item = numeric(1),
  top_item = numeric(1)
)

for (i in seq_len(nrow(bernoulli_data))) {
  r <- rankings[bernoulli_data$assessor[[i]], ]
  items_to_compare <- c(
    items[bernoulli_data$comp[[i]], 1],
    items[bernoulli_data$comp[[i]], 2]
  )

  a <- if (runif(1) > theta) {
    order(r[items_to_compare], decreasing = TRUE)
  } else {
    order(r[items_to_compare], decreasing = FALSE)
  }
  bernoulli_data[i, c("bottom_item", "top_item")] <- items_to_compare[a]
}


bernoulli_data$comp <- NULL
usethis::use_data(bernoulli_data, overwrite = TRUE)
