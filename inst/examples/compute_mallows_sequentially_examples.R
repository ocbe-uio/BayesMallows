# Observe one ranking at each of 12 timepoints
data <- lapply(seq_len(nrow(potato_visual)), function(i) {
  setup_rank_data(potato_visual[i, ])
})

