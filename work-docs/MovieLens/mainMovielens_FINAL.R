

library(BayesMallows)

compute_tc <- FALSE

if(compute_tc){
  load("./work-docs/movielens/movielensData_FINAL.RData")

  # prepare the data
  movie_tc <- generate_transitive_closure(ranking)
  movie_init_rank <- generate_initial_ranking(movie_tc)
  constraints <- generate_constraints(movie_tc, n_items = ncol(movie_init_rank))

  save(movie_tc, movie_init_rank, constraints, file = "./work-docs/movielens/movie_inits.RData")
} else {
  load("./work-docs/movielens/movie_inits.RData")
}

# compute the asymptotics partition function
alpha_vector <- c(seq(.01, 10, len = 50), seq(10.4, 40, len = 50))
degree <- 10
metric <- 'footrule'
n_iterations <- 500
K <- 80
n <- ncol(movie_init_rank)

estimate <- estimate_partition_function(method = "asymptotic",
                                        alpha_vector = alpha_vector,
                                        n_items = n, metric = metric,
                                        n_iterations = n_iterations,
                                        K = K, degree = degree)



# mcmc settings
n_clusters <- 1:20
nmc <- 1e3
burnin <- nmc-1
rho_thin <- nmc-1
aug_thin <- nmc-1
clus_ass_thin <- nmc - 1
L <- 20
sigmaAlpha <- .05
alphaJump <- 10
lambda <- .1

# mcmc
resList <- NULL
for(k in n_clusters){
  fitMallows <- compute_mallows(rankings = movie_init_rank, preferences = movie_tc,
                                nmc = nmc, n_clusters = k, save_augmented_data = FALSE,
                                rho_thinning = rho_thin, include_wcd = TRUE, logz_estimate = estimate,
                                leap_size = L, alpha_prop_sd = sigmaAlpha, alpha_jump = alphaJump,
                                lambda = lambda, aug_thinning = aug_thin, cluster_assignment_thinning = clus_ass_thin,
                                verbose = TRUE,
                                validate_rankings = FALSE, constraints = constraints,
                                skip_postprocessing = TRUE)

  print(paste('Finished cluster ',k,sep=''))
  save(fitMallows, file=paste('./work-docs/movielens/BayesMallows_nmc=',nmc,'_cluster',k,'.RData',sep=''))
  gc()
  #resList <- c(resList, list(fitMallows))
}



# pdf(file = paste('../../data/valeriv/results/Movielens/BayesMallows_nmc=',nmc,'_elbow.pdf',sep=''), width = 13, height = 7)
# plot_elbow(resList, burnin = burnin)
# dev.off()

  # this would be the correct way, but then the object becomes too big
  # models <- map(n_clusters, ~ compute_mallows(rankings = movie_init_rank, preferences = movie_tc,
  #                                             nmc = nmc, n_clusters = .x, include_wcd = TRUE,
  #                                             aug_diag_thinning = Rthin, thinning = thin))
  # save(models, file=paste('./Movielens/BayesMallows_nmc=',nmc,'_',as.character(max(n_clusters)),'clusters.RData',sep=''))
  # plot_elbow(models, burnin = burnin)


