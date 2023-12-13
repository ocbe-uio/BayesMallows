values <- read.delim("https://oeis.org/A126065/b126065.txt", sep = " ",
                     col.names = paste0("V", 1:2),
                     colClasses = c("NULL", "numeric"))$V2
n_items <- seq(from = 2, to = 60, by = 1)
ulam_cardinalities <- list(list(distance = 0, value = 1))
for(n in n_items) {
  distances <- seq(from = 0, to = n - 1, by = 1)
  if(length(distances) > length(values)) stop("error!")
  ulam_cardinalities[[n]] <- list(
    distance = distances,
    value = values[seq_along(distances)]
  )
  values <- values[-seq_along(distances)]
}

usethis::use_data(ulam_cardinalities, overwrite = TRUE)
