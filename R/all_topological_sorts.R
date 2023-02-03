# Translation to R of C++ and Python code found here
# https://www.geeksforgeeks.org/all-topological-sorts-of-a-directed-acyclic-graph/
all_topological_sorts <- function(graph, n_items, env, path = integer(),
                                  discovered = rep(FALSE, n_items)) {
  flag <- FALSE

  for (i in seq_len(n_items)) {
    if (attr(graph, "indegree")[[i]] == 0 && !discovered[[i]]) {
      attr(graph, "indegree")[graph[[i]]] <- attr(graph, "indegree")[graph[[i]]] - 1

      path <- c(path, i)
      discovered[[i]] <- TRUE
      all_topological_sorts(graph, n_items, env, path, discovered)

      attr(graph, "indegree")[graph[[i]]] <- attr(graph, "indegree")[graph[[i]]] + 1
      path <- path[-length(path)]
      discovered[[i]] <- FALSE

      flag <- TRUE
    }
  }
  if (length(path) == n_items) {
    assign("x", c(get("x", envir = env), list(path)), envir = env)
    assign("num", get("num", envir = env) + 1, envir = env)
  }
}
