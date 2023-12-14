values <- read.delim("https://oeis.org/A175929/b175929.txt",
  sep = " ",
  col.names = paste0("V", 1:2),
  colClasses = c("NULL", "numeric")
)$V2
n_items <- seq(from = 1, to = 20, by = 1)
values <- na.omit(values)
spearman_cardinalities <- list()
for (n in n_items) {
  distances <- seq(from = 0, to = 2 * choose(n + 1, 3), by = 2)
  if (length(distances) > length(values)) stop("error!")
  spearman_cardinalities[[n]] <- list(
    distance = distances,
    value = values[seq_along(distances)]
  )
  values <- values[-seq_along(distances)]
}

usethis::use_data(spearman_cardinalities, overwrite = TRUE)
