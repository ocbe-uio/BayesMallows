library(BayesMallows)
library(tidyverse)
library(parallel)

#load our data
data <- read.csv("work-docs/data.csv")

#we reshape it to fit BayesMallows package
flavours <- c(
  "Banan", "Caramel", "Chocolate Toffee", "Cocos", "Daim",
  "Eclairs", "Fransk nougat", "Golden Toffee", "Japp", "Lakris",
  "Marsipan", "Nougatcrisp", "NC8tti"
)
data_num <- apply(data[5:17], 2, function(x) match(x, flavours)) #replace characters with numbers
data_num_ord <- create_ordering(data_num) #replace ranking with ordering
colnames(data_num_ord) <- flavours #colnames and rownames just to remember
rownames(data_num_ord) <- paste0("Assessor", 1:dim(data_num)[1]) #colnames and rownames just to remember

data_BM <- setup_rank_data(rankings = data_num_ord)

#step 1 - observe and select burning
cl <- makeCluster(detectCores())
bmm <- compute_mallows_mixtures(
  n_clusters = c(1, 2, 3, 4),
  data = data_BM,
  compute_options = set_compute_options(nmc = 10000, include_wcd = FALSE),
  cl = cl)
stopCluster(cl)
assess_convergence(bmm)

selected_burnin <- 5000L
thinning <- 10L

#step 2 - observe and select the number of clusters
cl <- makeCluster(detectCores())
bmm <- compute_mallows_mixtures(
  n_clusters = 1:4,
  data = data_BM,
  compute_options =
    set_compute_options(nmc = 11000, burnin = selected_burnin,
                        rho_thinning = thinning, include_wcd = TRUE),
  cl = cl)
stopCluster(cl)
plot_elbow(bmm)


selected_n_clusters <- 3

#step 3 - compute clustering and see the results
bmm <- compute_mallows(
  data = data_BM,
  model_options = set_model_options(n_cluster = selected_n_clusters),
  compute_options = set_compute_options(
    nmc = 11000, burnin = selected_burnin,
    clus_thinning = thinning, rho_thinning = thinning, include_wcd = TRUE) #save_ind_clus = TRUE
)

assign_cluster(bmm)
plot(bmm, type = "cluster_probs")
