# running script to assess results

require(BayesMallows)
require(PerMallows)
require(dplyr)
require(tibble)
require(Rcpp)
require(ggplot2)
require(plotrix)
require(purrr)
require(crayon)
require(utf8)
require(fields)
require(tidyr)
require(dplyr)
require(gridExtra)
require(nimble)
######################
## Source Functions
######################

source("get_mallows_loglik.R")
source("metropolis_hastings_rho.R")
source("leap_and_shift_probs.R")
source("metropolis_hastings_alpha.R")
source("smc_mallows_new_users_complete.R")
source("post_processing_functions.R")

#########################
# Generate Dataset
#########################
set.seed(994)

load("sushi_rankings.rda")
data = sushi_rankings[1:100,]

# General
n_items <- dim(sushi_rankings)[2]  # Number of items
leap_size = floor(n_items/5)
metric = "footrule"

# Generate estimate of Z_n(alpha)
alpha_vector <- seq(from = 0, to = 15, by = 0.1)
iter = 1e4
degree <- 10

# Estimate the logarithm of the partition function of the Mallows rank model using the estimate partition function
logz_estimate <- estimate_partition_function(method = "importance_sampling",
                                             alpha_vector = alpha_vector,
                                             n_items = n_items, metric = metric,
                                             nmc = iter, degree = degree)


######################################
# BayesMallows Analysis (MCMC)
######################################
nmc = 10000
burnin=5000
model_fit <- compute_mallows(rankings = data, nmc = nmc, metric = metric, leap_size =leap_size, 
                             alpha_prop_sd = 0.15, logz_estimate = logz_estimate)

model_fit$burnin = burnin
assess_convergence(model_fit, parameter = "alpha")
assess_convergence(model_fit, parameter = "rho", items = 1:5)

plot(model_fit, parameter = "alpha")
plot(model_fit, parameter = "rho")

alpha_samples_table = data.frame(iteration = 1:nmc , value = model_fit$alpha$value)
alpha_samples_table = alpha_samples_table[(burnin+1):nmc,]

plot_posterior_alpha <- ggplot2::ggplot(alpha_samples_table, ggplot2::aes(x = alpha_samples_table$value)) +
  ggplot2::geom_density() +
  ggplot2::xlab(expression(alpha)) +
  scale_x_continuous(limits = c(1.3, 2.1)) +
  ggplot2::ylab("Posterior density") +
  ggplot2::ggtitle(label = "BayesMallows") + 
  theme(plot.title = element_text(hjust = 0.5)) 
print(plot_posterior_alpha)

# We can also compute the CP consensus posterior ranking
compute_consensus(model_fit, type = "CP")
compute_consensus(model_fit, type = "MAP")

# And we can compute the posterior intervals:
compute_posterior_intervals(model_fit, parameter = "alpha")
compute_posterior_intervals(model_fit, parameter = "rho")

# from observing the plots, this looks like the estimated parameters of the Mallows Model
rho_0 = c(4,5,2,6,8,3,9,1,7,10)
alpha_0 = 1.7

# heatplot - there is no burnin!
mcmc_rho_matrix = matrix(model_fit$rho$value, ncol = n_items, nrow = nmc, byrow=TRUE)
mcmc_heatmat_rho = heatMat(mcmcOutput = mcmc_rho_matrix, burnin = burnin, t_rank = rho_0)
mcmc_heatplot_full = heatPlot_fixed(mat = mcmc_heatmat_rho, t_rank = rho_0)


###################################################################
# SMC
###################################################################
mcmc_times = 5
num_new_obs = 10
Time = dim(data)[1]/num_new_obs
N = 100

test =  smc_mallows_new_users_complete(R_obs = data, n_items = n_items, metric = metric,
                                       leap_size = leap_size, N = N, Time = Time, 
                                       logz_estimate = logz_estimate, mcmc_kernel_app = mcmc_times, 
                                       num_new_obs = num_new_obs)

###############################
# Analysis
###############################

# posterior confidence intervals for rho
compute_posterior_intervals_rho(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0)
colnames(sushi_rankings)
# MAP AND CP consensus ranking estimates
rho_cp = compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "CP")
rho_cp
rho_map = compute_rho_consensus(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, C = 1, type = "MAP")
rho_map


# posterior for alpha
plot_alpha_posterior(output = test$alpha_samples[,Time+1], nmc = N, burnin = 0)

alpha_samples_table = data.frame(iteration = 1:N , value = test$alpha_samples[,Time+1])
plot_posterior_alpha <- ggplot2::ggplot(alpha_samples_table, ggplot2::aes(x = alpha_samples_table$value)) +
  ggplot2::geom_density() +
  ggplot2::xlab(expression(alpha)) +
  scale_x_continuous(limits = c(1.3, 2.1)) +
  ggplot2::ylab("Posterior density") +
  ggplot2::ggtitle(label = "Implemented SMC scheme") + 
  theme(plot.title = element_text(hjust = 0.5)) 
print(plot_posterior_alpha)

# posterior confidence intervals
alpha_posterior_intervals = compute_posterior_intervals_alpha(output = test$alpha_samples[,Time+1], nmc = N, burnin = 0)

# heatplot for rho
plot_rho_heatplot(output = test$rho_samples[,,Time+1], nmc = N, burnin = 0, n_items = n_items, rho_0 = rho_0)
