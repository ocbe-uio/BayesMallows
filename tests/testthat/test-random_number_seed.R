context("Testing that the random numbers are equal between platforms")

bmm <- compute_mallows_mixtures(
  n_clusters = c(1, 4),
  rankings = sushi_rankings,
  compute_options = set_compute_options(nmc = 50, include_wcd = FALSE),
  seed = 432
)

expect_equal(
  bmm[[2]]$cluster_probs$value[20:26],
  c(
    0.103780617148549, 0.376634968419074, 0.133399257941488, 0.408615632333589,
    0.0813501413058496, 0.409409209490965, 0.110568949260592
  )
)


m <- compute_mallows(sushi_rankings, n_clusters = 5, seed = 123,
                     compute_options = set_compute_options(nmc = 20))

expect_equal(
  m$cluster_probs$value[30:35],
  c(
    0.735384849956839, 0.113779597344025, 0.0230786398725923, 0.022571383303261,
    0.052670618132848, 0.787899761347273
  )
)
