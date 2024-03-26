devtools::load_all()

preferences = subset(beach_preferences, bottom_item < 10 & top_item < 10)
random = TRUE
shuffle_unranked = TRUE
random_limit = 100

prefs <- splitpref(preferences)
n_items <- 10


mat <- prefs[[1]]

system.time({
  all_topological_sorts_cpp(mat, n_items)
})

system.time({
  graph <- list()
  for (i in seq_len(n_items)) {
    graph[[i]] <- unique(mat[mat[, "top_item"] == i, "bottom_item"])
  }
  indegree_init <- rep(0, n_items)
  indegree <- table(unlist(graph))
  indegree_init[as.integer(names(indegree))] <- indegree
  attr(graph, "indegree") <- indegree_init

  e1 <- new.env()
  assign("x", list(), envir = e1)
  assign("num", 0L, envir = e1)
  all_topological_sorts(graph, n_items, e1)

  get("x", envir = e1)[[sample(get("num", envir = e1), 1)]]

})

do.call(rbind, get("x", envir = e1))
