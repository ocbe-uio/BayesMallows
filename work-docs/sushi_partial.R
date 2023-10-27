rm(list = ls())
devtools::load_all()
data_partial <- sushi_rankings
data_partial[data_partial > 5] <- NA

data_batch1 <- data_partial[1:300, ]
model <- compute_mallows(rankings = data_batch1)
assess_convergence(model)
model$burnin <- 300

compute_consensus(model, type = "MAP")

new_rankings <- data_partial[301:600, ]
model2 <- update_mallows(model, new_rankings = new_rankings,
                         n_particles = 1000, type = "partial")

compute_consensus(model2)
compute_consensus(model2, type = "MAP")

new_rankings <- sushi_rankings[601:900, ]
model3 <- update_mallows.SMCMallows(model = model2, new_rankings = new_rankings)

compute_consensus(model3)
compute_consensus(model3, type = "MAP")
