context('SMC-mallows new users partial rankings')

#########################
# Generate Dataset
#########################
set.seed(994)

# General
n_items <- dim(sushi_rankings)[2]  # Number of items
rho_0 <- seq(from = 1, to = n_items, by = 1) # 'true' consensus ranking
alpha_0 <- 2  # fixed/ 'true' scale parameter
leap_size = floor(n_items/5)
metric = "footrule"

# Generate estimate of Z_n(alpha)
alpha_vector <- seq(from = 0, to = 20, by = 0.1)
iter = 1e4
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(method = "importance_sampling",
                                             alpha_vector = alpha_vector,
                                             n_items = n_items, metric = metric,
                                             nmc = iter, degree = degree)

# Make this information partially observed over time
samples = sushi_rankings[1:100,]
samples[samples > 5] <- NA


#######################################
# Bayesmallows MCMC Results
#######################################
nmc = 2000
bayesmallows_mcmc <- compute_mallows(samples, nmc = nmc, leap_size =leap_size, metric = metric, alpha_prop_sd = 0.15)
bayesmallows_mcmc$burnin = 1000

# choice items to see in trace plot
items <- sample(1:n_items, 5, replace=F)
items = sort(items)

# plot trace - rho
# p <- assess_convergence(bayesmallows_mcmc, parameter = "rho", items = items)  +
#   ggtitle("Trace plot for MCMC algorithm in BayesMallows R package") +
#   theme(plot.title = element_text(hjust = 0.5))
# p

# heatplot
bayesmallows_rho_matrix = matrix(bayesmallows_mcmc$rho$value, ncol = n_items, nrow = nmc, byrow=TRUE)
bayesmallows_heatmat_rho = BayesMallows:::heatMat(mcmcOutput = bayesmallows_rho_matrix, burnin = bayesmallows_mcmc$burnin, t_rank = rho_0)
# bayesmallows_heatplot = BayesMallows:::heatPlot_fixed(mat = bayesmallows_heatmat_rho,t_rank = rho_0)


# posterior intervals for rho
compute_posterior_intervals(bayesmallows_mcmc, parameter = "rho")

compute_consensus(model_fit = bayesmallows_mcmc, type = "CP", burnin = bayesmallows_mcmc$burnin)
compute_consensus(model_fit = bayesmallows_mcmc, type = "MAP", burnin = bayesmallows_mcmc$burnin)

#plot trace - alpha
# assess_convergence(bayesmallows_mcmc, parameter = "alpha")+
#   ggtitle("Trace plot for MCMC algorithm in BayesMallows R package") +
#   theme(plot.title = element_text(hjust = 0.5))
#ggsave("bayesmallows_alpha_traceplot.png", height = 6, width = 7)

# bayesmallows_alpha_posterior_plot = plot(bayesmallows_mcmc, parameter = "alpha", burnin = bayesmallows_mcmc$burnin) +
#   ggtitle("BayesMallows R package") + theme(plot.title = element_text(hjust = 0.5))
# bayesmallows_alpha_posterior_plot
#ggsave("bayesmallows_alpha_posterior_plot.png", height = 6, width = 7)

# determine posterior estimates - alpha
bayesmallows_alpha_posterior_intervals = compute_posterior_intervals(bayesmallows_mcmc, parameter = "alpha")
bayesmallows_alpha_posterior_intervals




###############################
# SMC Analysis (alpha unknown)
###############################

N = 100
mcmc_times = 5
num_new_obs = 5
Time = dim(samples)[1]/num_new_obs

alpha_prop_sd = 0.5
lambda = 0.15
alpha_max = 1e6
aug_method = "random"
#aug_method = "pseudolikelihood"

test =  smc_mallows_new_users_partial(R_obs = samples, n_items = n_items, metric = metric,
                                             leap_size = leap_size, N = N, Time = Time,
                                             logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times,
                                             num_new_obs = num_new_obs, alpha_prop_sd = alpha_prop_sd,
                                             lambda = lambda, alpha_max = alpha_max, aug_method = aug_method)


# heatplot for rho

# plot_rho_heatplot(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, n_items = n_items, rho_0 = rho_0)

# posterior confidence intervals for rho

BayesMallows:::compute_posterior_intervals_rho(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0)

# MAP AND CP consensus ranking estimates
rho_cp = BayesMallows:::compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "CP")
rho_cp
rho_map = BayesMallows:::compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "MAP")
rho_map


# posterior for alpha
# plot_alpha_posterior(output = test$alpha_samples[,Time+1], nmc = N, burnin = 0)

# posterior confidence intervals
alpha_posterior_intervals = BayesMallows:::compute_posterior_intervals_alpha(output = test$alpha_samples[,Time+1], nmc = N, burnin = 0)




###############################
# SMC Analysis (alpha fixed)
###############################

alpha_0 = 2

test =  smc_mallows_new_users_partial_alpha_fixed(R_obs = samples, n_items = n_items, metric = metric,
                                      leap_size = leap_size, N = N, Time = Time,
                                      logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times,
                                      num_new_obs = num_new_obs, aug_method = aug_method, alpha = alpha_0)

# heatplot for rho
# plot_rho_heatplot_partial(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, n_items = n_items, rho_0 = rho_0)

# posterior confidence intervals for rho
compute_posterior_intervals_rho(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0)

# MAP AND CP consensus ranking estimates
rho_cp = compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "CP")
rho_cp
rho_map = compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "MAP")
rho_map


