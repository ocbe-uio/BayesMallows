library(dplyr)
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4
)

pair_comp_tc <- generate_transitive_closure(pair_comp)

initial_ranking <- generate_initial_ranking(pair_comp_tc)

# set.seed(1)
# res <- compute_mallows(initial_ranking, pair_comp_tc, nmc = 10)



rankings = initial_ranking
preferences = pair_comp_tc
metric = "footrule"
n_clusters = 1L
nmc = 2L
leap_size = NULL
rho_init = NULL
thinning = 1L
alpha_prop_sd = 0.1
alpha_init = 1
alpha_jump = 1L
lambda = 0.001
psi = 10L
include_wcd = (n_clusters > 1)
save_augmented_data = FALSE
is_fit = NULL
devtools::load_all()



stopifnot(!is.null(rankings) || !is.null(preferences))

if(!is.null(rho_init)) {
  stopifnot(validate_permutation(rho_init) && (sum(is.na(rho_init)) == 0))
}

# Deal with pairwise comparisons. Generate rankings compatible with them.
if(!is.null(preferences)){

  if(!("BayesMallowsTC" %in% class(preferences))){
    preferences <- generate_transitive_closure(preferences)
  }

  if(is.null(rankings)){
    rankings <- generate_initial_ranking(as.matrix(preferences))
  }

  linear_ordering <- split(preferences, preferences$assessor)
  linear_ordering <- purrr::map(linear_ordering, function(x) {
    create_linear_ordering(as.matrix(x[,2:3]), partial = TRUE)
  })

} else {
  linear_ordering <- list()
}

# Check that all rows of rankings are proper permutations
if(!all(apply(rankings, 1, validate_permutation))){
  stop("Not valid permutation.")
}

# Check that we do not jump over all alphas
stopifnot(alpha_jump < nmc)

# Check that we do not jump over all rhos
stopifnot(thinning < nmc)

# Find the number of items
n_items <- ncol(rankings)

# Set leap_size if it is not alredy set.
if(is.null(leap_size)) leap_size <- floor(n_items / 5)

# Extract the right sequence of cardinalities, if relevant
if(!is.null(is_fit)){
  cardinalities <- NULL
  message("Using user-provided importance sampling estimate of partition function.")
} else if(metric %in% c("footrule", "spearman")){
  # Extract the relevant rows from partition_function_data
  # Note that we need to evaluate the right-hand side, in particular metric,
  # to avoid confusion with columns of the tibble
  relevant_params <- dplyr::filter(partition_function_data,
                                   .data$n_items == !!n_items,
                                   .data$metric == !!metric
  )

  type <- dplyr::pull(relevant_params, type)

  if(type == "cardinalities") {
    cardinalities <- unlist(relevant_params$values)
    is_fit <- NULL
  } else if(type == "importance_sampling"){
    cardinalities <- NULL
    is_fit <- unlist(relevant_params$values)
  } else {
    stop("Precomputed partition function not available yet. Consider computing one
           with the function estimate_partition_function(), and provide it
           in the is_fit argument to compute_mallows().")
  }

} else if (metric %in% c("cayley", "hamming", "kendall")) {
  cardinalities <- NULL
  is_fit <- NULL
} else {
  stop(paste("Unknown metric", metric))
}


element <- 1
BayesMallows:::find_pairwise_limits(0, 6, element, linear_ordering[[1]],
                                    possible_rankings = initial_ranking[1, ]) + c(1, -1)
BayesMallows:::find_pairwise_limits_old(0, 6, element, as.matrix(preferences),
                                        possible_rankings = initial_ranking[1, ], 0)





# to extract one sample at a time. armadillo is column major, just like rankings
fit <- run_mcmc(rankings = t(rankings),
                nmc = nmc,
                linear_ordering = linear_ordering,
                cardinalities = cardinalities,
                is_fit = is_fit,
                rho_init = rho_init,
                metric = metric,
                n_clusters = n_clusters,
                include_wcd = include_wcd,
                lambda = lambda,
                leap_size = leap_size,
                alpha_prop_sd = alpha_prop_sd,
                alpha_init = alpha_init,
                alpha_jump = alpha_jump,
                thinning = thinning,
                save_augmented_data = save_augmented_data
)
