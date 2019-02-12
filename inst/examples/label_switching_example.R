\dontrun{
  # This example shows how to assess if label switching happens in BayesMallows

  library(BayesMallows)
  # We start by creating a directory in which csv files with individual
  # cluster probabilities should be saved in each step of the MCMC algorithm
  dir.create("./test_label_switch")
  # Next, we go into this directory
  setwd("./test_label_switch/")
  # For comparison, we run compute_mallows with and without saving the cluster
  # probabilities The purpose of this is to assess the time it takes to save
  # the cluster probabilites
  system.time(m <- compute_mallows(rankings = sushi_rankings,
                                   n_clusters = 6, nmc = 2000, save_clus = TRUE,
                                   save_ind_clus = FALSE))
  # With this options, compute_mallows will save cluster_probs2.csv,
  # cluster_probs3.csv, ..., cluster_probs[nmc].csv.
  system.time(m <- compute_mallows(rankings = sushi_rankings, n_clusters = 6,
                                   nmc = 2000, save_clus = TRUE,
                                   save_ind_clus = TRUE))

  # Next, we check convergence of alpha
  assess_convergence(m)

  # We set the burnin to 1000
  burnin <- 1000

  # Find all files that were saved. Note that the first file saved is cluster_probs2.csv
  cluster_files <- list.files(pattern = "cluster\\_probs[[:digit:]]+\\.csv")

  # Check the size of the files that were saved.
  paste(sum(do.call(file.size, list(cluster_files))) * 1e-6, "MB")

  # Find the iteration each file corresponds to, by extracting its number
  library(stringr)
  iteration_number <- as.integer(str_extract(cluster_files, "[:digit:]+"))
  # Remove all files before burnin
  file.remove(cluster_files[iteration_number <= burnin])
  # Update the vector of files, after the deletion
  cluster_files <- list.files(pattern = "cluster\\_probs[[:digit:]]+\\.csv")
  # Create 3d array, with dimensions (iterations, assessors, clusters)
  prob_array <- array(dim = c(length(cluster_files), m$n_assessors, m$n_clusters))
  # Read each file, adding to the right element of the array
  library(readr)
  for(i in seq_along(cluster_files)){
    prob_array[i, , ] <- as.matrix(
      read_delim(cluster_files[[i]], delim = ",",
                 col_names = FALSE, col_types = paste(rep("d", m$n_clusters),
                                                      collapse = "")))
  }

  library(dplyr)
  library(tidyr)
  # Create an tnteger array of latent allocations, as this is required by label.switching
  z <- m$cluster_assignment %>%
    filter(iteration > burnin) %>%
    mutate(value = as.integer(str_extract(value, "[:digit:]+"))) %>%
    spread(key = assessor, value = value, sep = "_") %>%
    select(-iteration) %>%
    as.matrix()

  # Now apply Stephen's algorithm
  library(label.switching)
  ls <- label.switching("STEPHENS", z = z, K = m$n_clusters, p = prob_array)

  # Check the proportion of cluster assignments that were switched
  mean(apply(ls$permutations$STEPHENS, 1, function(x) !all.equal(x, seq(1, m$n_clusters))))

  # Remove the rest of the csv files
  file.remove(cluster_files)
  # Move up one directory
  setwd("..")
  # Remove the directory in which the csv files were saved
  file.remove("./test_label_switch/")
}
