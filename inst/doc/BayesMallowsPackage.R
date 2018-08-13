## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
set.seed(232312)

## ------------------------------------------------------------------------
library(BayesMallows)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(potato_visual, caption = "Example dataset `potato_visual`.")

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual)

## ---- cache=TRUE---------------------------------------------------------
str(model_fit)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit, type = "rho", items = 1:5)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit, type = "rho", items = c("P16", "P17", "P18", "P19", "P20"))

## ---- cache=TRUE---------------------------------------------------------
plot(model_fit, burnin = 2000)

## ---- fig.show='hold', cache=TRUE----------------------------------------
model_fit_big <- compute_mallows(potato_visual, nmc = 1e4 + 2000)

## ---- cache=TRUE---------------------------------------------------------
plot(model_fit_big, burnin = 2000)

## ---- cache=TRUE---------------------------------------------------------
rm(model_fit_big)

## ---- fig.width=6, fig.height=6, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- cache=TRUE---------------------------------------------------------
test_run <- compute_mallows(potato_visual, nmc = 10000, alpha_jump = 10)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(test_run, type = "alpha")

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(test_run, type = "rho", items = 1:5)

## ---- cache=TRUE---------------------------------------------------------
rm(test_run)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual, nmc = 25000 + 2000, alpha_jump = 10)

## ---- cache=TRUE---------------------------------------------------------
plot(model_fit, burnin = 2000)

## ---- fig.width=6, fig.height=6, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 2)

## ---- fig.width=6, fig.height=6, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- cache=TRUE---------------------------------------------------------
t1 <- Sys.time()
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 1)
t2 <- Sys.time()
(no_thinning <- t2 - t1)
t3 <- Sys.time()
model_fit <- compute_mallows(potato_visual, nmc = 50000 + 2000, 
                             alpha_jump = 10, thinning = 10)
t4 <- Sys.time()
(thinning <- t4 - t3)


## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual, metric = "kendall",
                             nmc = 25000 + 2000, alpha_jump = 10)

## ---- fig.width=6, fig.height=6, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual, metric = "spearman",
                             nmc = 25000 + 2000, alpha_jump = 10)

## ---- fig.width=6, fig.height=6, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- error=TRUE, cache=TRUE---------------------------------------------
potato_modified <- potato_visual
potato_modified[1, 1:2] <- 1

model_fit <- compute_mallows(potato_modified)

## ---- message=FALSE, cache=TRUE------------------------------------------
library(dplyr)
potato_top <- potato_visual * if_else(potato_visual > 5, NA_integer_, 1L)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(potato_top, caption = "Example dataset potato_top.")

## ------------------------------------------------------------------------
item_ranked <- apply(potato_top, 2, function(x) !all(is.na(x)))
potato_top <- potato_top[, item_ranked, drop = FALSE]

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(potato_top, caption = "Example dataset `potato_top`.")

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_top)

## ---- cache=TRUE---------------------------------------------------------
str(model_fit)

## ---- fig.width=6, fig.height=5, cache=TRUE------------------------------
assess_convergence(model_fit, type = "augmentation")

## ---- cache=TRUE---------------------------------------------------------
model_fit1 <- compute_mallows(potato_top, nmc = 10000)
model_fit2 <- compute_mallows(potato_top, nmc = 10000)

## ---- fig.show='hold', cache=TRUE----------------------------------------
assess_convergence(model_fit1, type = "augmentation", assessors = 3)
assess_convergence(model_fit2, type = "augmentation", assessors = 3)

## ---- cache=TRUE---------------------------------------------------------
rm(model_fit, model_fit2)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit1, type = "alpha")

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit1, type = "rho", items = 1:8)

## ---- cache=TRUE---------------------------------------------------------
rm(model_fit1)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_top, nmc = 1e4 + 2e4)

## ---- fig.width=6, cache=TRUE--------------------------------------------
plot(model_fit, burnin = 1e4)

## ---- fig.width=6, fig.height=5, cache=TRUE------------------------------
plot(model_fit, burnin = 1e4, type = "rho", items = 1:8)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_top, nmc = 1e4 + 2e4, metric = "spearman")

## ---- fig.width=6, fig.height=5, cache=TRUE------------------------------
plot(model_fit, burnin = 1e4, type = "rho", items = 1:8)

## ---- cache=TRUE---------------------------------------------------------
missing_indicator <- if_else(
  runif(nrow(potato_visual) * ncol(potato_visual)) < 0.1,
                            NA_real_, 1)
potato_missing <- potato_visual * missing_indicator

## ---- echo=FALSE, results='asis', cache=TRUE-----------------------------
knitr::kable(potato_missing, caption = "Example dataset `potato_missing`.")

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_missing, nmc = 1e4)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit, type = "rho", items = 1:6)

## ---- fig.width=6, fig.height=5, cache=TRUE------------------------------
assess_convergence(model_fit, type = "augmentation")

## ---- cache=TRUE---------------------------------------------------------
apply(potato_missing, 1, function(x) sum(is.na(x)))

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(potato_visual, nmc = 1e4)

## ---- fig.width=6, fig.height=5, cache=TRUE------------------------------
plot(model_fit, burnin = 2000, type = "rho", items = 1:20)

## ---- message=FALSE, results='asis', cache=TRUE--------------------------
pair_comp <- tribble(
  ~assessor, ~bottom_item, ~top_item,
  1, 1, 2,
  1, 2, 5,
  1, 4, 5,
  2, 1, 2,
  2, 2, 3,
  2, 3, 4
)

knitr::kable(pair_comp)

## ---- cache=TRUE---------------------------------------------------------
pair_comp_tc <- generate_transitive_closure(pair_comp)

## ---- results='asis', cache=TRUE-----------------------------------------
knitr::kable(pair_comp_tc)

## ---- cache=TRUE---------------------------------------------------------
class(pair_comp_tc)

## ---- cache=TRUE---------------------------------------------------------
initial_ranking <- generate_initial_ranking(pair_comp_tc)

## ---- results='asis', cache=TRUE-----------------------------------------
knitr::kable(initial_ranking, row.names = TRUE)

## ---- cache=TRUE---------------------------------------------------------
model_fit <- compute_mallows(R = initial_ranking, P = pair_comp_tc)

## ---- cache=TRUE---------------------------------------------------------
str(model_fit)

## ---- fig.width=6, cache=TRUE--------------------------------------------
assess_convergence(model_fit, type = "augmentation")

## ---- results='asis'-----------------------------------------------------
knitr::kable(head(beach_preferences, 6), caption = "Example dataset `beach_preferences`")

## ---- cache=TRUE---------------------------------------------------------
beach_tc <- generate_transitive_closure(beach_preferences)

## ---- cache=TRUE---------------------------------------------------------
str(beach_preferences)
str(beach_tc)

## ---- cache=TRUE---------------------------------------------------------
beach_init_rank <- generate_initial_ranking(beach_tc)

## ---- cache=TRUE, results='asis'-----------------------------------------
knitr::kable(head(beach_init_rank, 6))

## ---- error=TRUE---------------------------------------------------------
test_run <- compute_mallows(R = beach_init_rank, P = beach_tc, nmc = 2)

