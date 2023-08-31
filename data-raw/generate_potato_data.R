# This script generates the potato data based on potato_data.RData
# The RData file comes from the original data in Oystein Sorensen's PhD thesis
load("./data-raw/potato_data.RData")

potato_visual <- visualRo
colnames(potato_visual) <- paste0("P", seq_len(ncol(potato_visual)))
rownames(potato_visual) <- paste0("A", seq_len(nrow(potato_visual)))

potato_weighing <- weighingRo
colnames(potato_weighing) <- paste0("P", seq_len(ncol(potato_weighing)))
rownames(potato_weighing) <- paste0("A", seq_len(nrow(potato_weighing)))

potato_true_ranking <- Rtrue
names(potato_true_ranking) <- paste0("P", seq_along(potato_true_ranking))


devtools::use_data(potato_visual, potato_weighing, potato_true_ranking,
  overwrite = TRUE
)
