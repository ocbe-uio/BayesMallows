# Testing with footrule
set.seed(200)
model_fit <- compute_mallows(setup_rank_data(potato_weighing),
  model = set_model_options(metric = "footrule"),
  compute_options = set_compute_options(nmc = 1000)
)
mean_alpha <- mean(model_fit$alpha$value[501:1000])

test_that(
  "alpha is in a decent range for footrule",
  {
    expect_true(mean_alpha > 10 && mean_alpha < 20)
  }
)





test_that(
  "rho is a rank vector",
  {
    expect_true({
      obj <- aggregate(list(n = seq_len(nrow(model_fit$rho))),
        by = model_fit$rho, FUN = length
      )
      nrow(obj[obj$n > 1, ]) == 0
    })
  }
)
