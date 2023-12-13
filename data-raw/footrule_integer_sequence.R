values <- read.delim("https://oeis.org/A062869/b062869.txt", sep = " ")$X1
n_items <- seq(from = 1, to = 50, by = 1)

footrule_cardinalities <- list()
for(n in n_items) {
  distances <- seq(from = 0, to = floor(n^2 / 2), by = 2)
  footrule_cardinalities[[n]] <- list(
    distance = distances,
    value = values[seq_along(distances)]
  )
  values <- values[-seq_along(distances)]
}

usethis::use_data(footrule_cardinalities, internal = TRUE, overwrite = TRUE)
