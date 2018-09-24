

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
n_clusters <- 20:1
nmc <- 1e5
burnin <- nmc-1
rho_thin <- nmc-1
aug_thin <- nmc-1
clus_ass_thin <- nmc - 1
L <- 20
sigmaAlpha <- .05
alphaJump <- 10
lambda <- .1

# mcmc
models <- compute_mallows_mixtures(n_clusters = n_clusters,
                                   rankings = movie_init_rank,
                                   preferences = movie_tc,
                                   nmc = nmc,
                                   save_augmented_data = FALSE,
                                   rho_thinning = rho_thin,
                                   logz_estimate = estimate,
                                   leap_size = L,
                                   alpha_prop_sd = sigmaAlpha,
                                   alpha_jump = alphaJump,
                                   lambda = lambda,
                                   aug_thinning = aug_thin,
                                   cluster_assignment_thinning = clus_ass_thin,
                                   verbose = TRUE,
                                   validate_rankings = FALSE,
                                   constraints = constraints,
                                   skip_postprocessing = TRUE)

# Tidy the WCD
# NB! This is just my way of using an internal function. We will decide on the user interface later.
models <- purrr::map(models, BayesMallows:::tidy_wcd)

plot_elbow(models, burnin = floor(nmc/2))
ggplot2::ggsave(file = "elbow_20_clusters.png")
# pdf(file = paste('../../data/valeriv/results/Movielens/BayesMallows_nmc=',nmc,'_elbow.pdf',sep=''), width = 13, height = 7)
# plot_elbow(resList, burnin = burnin)
# dev.off()

  # this would be the correct way, but then the object becomes too big
  # models <- map(n_clusters, ~ compute_mallows(rankings = movie_init_rank, preferences = movie_tc,
  #                                             nmc = nmc, n_clusters = .x, include_wcd = TRUE,
  #                                             aug_diag_thinning = Rthin, thinning = thin))
  # save(models, file=paste('./Movielens/BayesMallows_nmc=',nmc,'_',as.character(max(n_clusters)),'clusters.RData',sep=''))
  # plot_elbow(models, burnin = burnin)


